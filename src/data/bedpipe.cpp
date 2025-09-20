#include "gelex/data/bedpipe.h"

#include <expected>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader.h"
#include "gelex/error.h"
#include "gelex/logger.h"

namespace gelex
{

BedPipe::BedPipe(
    std::ifstream&& file_stream,
    std::unique_ptr<detail::FamLoader> fam_loader,
    std::unique_ptr<detail::BimLoader> bim_loader,
    Eigen::Index bytes_per_variant,
    std::filesystem::path bed_path)
    : file_stream_(std::move(file_stream)),
      fam_loader_(std::move(fam_loader)),
      bim_loader_(std::move(bim_loader)),
      bytes_per_variant_(bytes_per_variant),
      bed_path_(std::move(bed_path))
{
    // File is already validated and opened by create()
    file_stream_.seekg(3);  // Skip the 3-byte magic number

    // Initialize both sample maps with original FAM order
    const auto& ids = fam_loader_->sample_ids();
    Eigen::Index index = 0;
    for (const auto& id : ids)
    {
        raw_sample_map_[id] = index;
        load_sample_map_[id] = index++;
    }
}

BedPipe::~BedPipe()
{
    if (file_stream_.is_open())
    {
        file_stream_.close();
    }
}

auto BedPipe::create(const std::filesystem::path& prefix, bool iid_only)
    -> std::expected<BedPipe, Error>
{
    // Create file paths
    std::filesystem::path bed_path = prefix;
    bed_path.replace_extension(".bed");
    std::filesystem::path fam_path = prefix;
    fam_path.replace_extension(".fam");
    std::filesystem::path bim_path = prefix;
    bim_path.replace_extension(".bim");

    // Load fam and bim files
    auto fam_loader = detail::FamLoader::create(fam_path.string(), iid_only);
    if (!fam_loader)
    {
        return std::unexpected(fam_loader.error());
    }

    auto bim_loader = detail::BimLoader::create(bim_path.string());
    if (!bim_loader)
    {
        return std::unexpected(bim_loader.error());
    }

    // Open and validate BED file
    auto file
        = detail::openfile<std::ifstream>(bed_path, detail::file_type::binary);
    if (!file)
    {
        return std::unexpected(file.error());
    }

    if (auto validation = validate_bed_file(*file, bed_path); !validation)
    {
        return std::unexpected(validation.error());
    }

    const auto num_samples
        = static_cast<Eigen::Index>(fam_loader->sample_ids().size());
    const Eigen::Index bytes_per_variant
        = calculate_bytes_per_variant(num_samples);

    auto logger = gelex::logging::get();
    logger->info(
        "Loaded BED file with {} samples and {} variants",
        fam_loader->sample_ids().size(),
        bim_loader->snp_ids().size());

    return BedPipe(
        std::move(*file),
        std::make_unique<detail::FamLoader>(std::move(*fam_loader)),
        std::make_unique<detail::BimLoader>(std::move(*bim_loader)),
        bytes_per_variant,
        std::move(bed_path));
}

auto BedPipe::validate_bed_file(
    std::ifstream& file,
    const std::filesystem::path& path) -> std::expected<void, Error>
{
    // Check magic number: 0x6c, 0x1b, 0x01
    uint8_t magic[3];
    file.read(reinterpret_cast<char*>(magic), 3);

    if (!file.good() || file.gcount() != 3)
    {
        return std::unexpected(enrich_with_file_info(
            Error{
                ErrorCode::InvalidFile, "Failed to read BED file magic number"},
            path.string()));
    }

    if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01)
    {
        return std::unexpected(enrich_with_file_info(
            Error{ErrorCode::InvalidFile, "Invalid BED file magic number"},
            path.string()));
    }

    return {};
}

auto BedPipe::calculate_bytes_per_variant(Eigen::Index num_samples)
    -> Eigen::Index
{
    return static_cast<Eigen::Index>(
        (num_samples + 3) / 4);  // Round up to nearest multiple of 4 samples
}

Eigen::Index BedPipe::calculate_offset(Eigen::Index variant_index) const
{
    return 3 + (variant_index * bytes_per_variant_);
}

auto BedPipe::validate_variant_index(Eigen::Index variant_index) const
    -> std::expected<void, Error>
{
    if (variant_index >= num_variants())
    {
        return std::unexpected(
            Error{ErrorCode::InvalidRange, "Variant index out of range"});
    }
    return {};
}

Eigen::VectorXd BedPipe::reorder_genotypes(
    const Eigen::VectorXd& raw_genotypes,
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
{
    Eigen::VectorXd reordered(id_map.size());

    for (const auto& [sample_id, load_index] : id_map)
    {
        const Eigen::Index raw_index = raw_sample_map_.at(sample_id);
        reordered(load_index) = raw_genotypes(raw_index);
    }

    return reordered;
}

auto BedPipe::get_genotypes(Eigen::Index variant_index) const
    -> std::expected<Eigen::VectorXd, Error>
{
    if (auto validation = validate_variant_index(variant_index); !validation)
    {
        return std::unexpected(validation.error());
    }

    const Eigen::Index offset = calculate_offset(variant_index);
    file_stream_.seekg(offset);

    std::vector<uint8_t> buffer(bytes_per_variant_);
    file_stream_.read(
        reinterpret_cast<char*>(buffer.data()), bytes_per_variant_);

    if (!file_stream_.good() || file_stream_.gcount() != bytes_per_variant_)
    {
        return std::unexpected(enrich_with_file_info(
            Error{
                ErrorCode::FileIOError,
                "Failed to read genotype data from BED file"},
            bed_path_.string()));
    }

    const auto num_samples
        = static_cast<Eigen::Index>(fam_loader_->sample_ids().size());
    Eigen::VectorXd genotypes(num_samples);

    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        const Eigen::Index byte_index = i / 4;
        const Eigen::Index bit_offset = 2 * (i % 4);
        const uint8_t genotype_code = (buffer[byte_index] >> bit_offset) & 0x03;
        genotypes(i) = add_map_[genotype_code];
    }

    return genotypes;
}

auto BedPipe::get_genotype(
    Eigen::Index variant_index,
    Eigen::Index sample_index) const -> std::expected<double, Error>
{
    const auto num_samples
        = static_cast<Eigen::Index>(fam_loader_->sample_ids().size());

    if (auto validation = validate_variant_index(variant_index); !validation)
    {
        return std::unexpected(validation.error());
    }

    if (sample_index >= num_samples)
    {
        return std::unexpected(
            Error{ErrorCode::InvalidRange, "Sample index out of range"});
    }

    const Eigen::Index offset
        = calculate_offset(variant_index) + (sample_index / 4);
    file_stream_.seekg(offset);

    uint8_t byte;
    file_stream_.read(reinterpret_cast<char*>(&byte), 1);

    if (!file_stream_.good())
    {
        return std::unexpected(enrich_with_file_info(
            Error{
                ErrorCode::FileIOError,
                "Failed to read genotype data from BED file"},
            bed_path_.string()));
    }

    const Eigen::Index bit_offset = 2 * (sample_index % 4);
    const uint8_t genotype_code = (byte >> bit_offset) & 0x03;

    return add_map_[genotype_code];
}

auto BedPipe::get_sample_genotypes(Eigen::Index sample_index) const
    -> std::expected<Eigen::VectorXd, Error>
{
    if (sample_index >= raw_sample_ids().size())
    {
        return std::unexpected(
            Error{ErrorCode::InvalidRange, "Sample index out of range"});
    }

    Eigen::VectorXd genotypes(num_variants());

    for (Eigen::Index i = 0; i < num_variants(); ++i)
    {
        auto result = get_genotype(i, sample_index);
        if (!result)
        {
            return std::unexpected(result.error());
        }
        genotypes(i) = *result;
    }

    return genotypes;
}

auto BedPipe::read_variants_bulk(
    Eigen::Index start_variant,
    Eigen::Index end_variant,
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
    -> std::expected<Eigen::MatrixXd, Error>
{
    const auto num_samples
        = static_cast<Eigen::Index>(fam_loader_->sample_ids().size());
    const Eigen::Index num_variants_in_chunk = end_variant - start_variant;

    if (num_variants_in_chunk <= 0)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidRange,
                "Invalid variant range for bulk reading"});
    }

    // Calculate the total bytes to read
    const Eigen::Index start_offset = calculate_offset(start_variant);
    const Eigen::Index total_bytes = num_variants_in_chunk * bytes_per_variant_;

    // Seek to the start position
    file_stream_.seekg(start_offset);

    // Read all variants in one go
    std::vector<uint8_t> buffer(total_bytes);
    file_stream_.read(reinterpret_cast<char*>(buffer.data()), total_bytes);

    if (!file_stream_.good() || file_stream_.gcount() != total_bytes)
    {
        return std::unexpected(enrich_with_file_info(
            Error{
                ErrorCode::FileIOError,
                "Failed to read genotype data from BED file in bulk"},
            bed_path_.string()));
    }

    // Parse all variants from the buffer (in raw sample order)
    Eigen::MatrixXd raw_genotypes(num_samples, num_variants_in_chunk);

    for (Eigen::Index variant_idx = 0; variant_idx < num_variants_in_chunk;
         ++variant_idx)
    {
        const uint8_t* variant_data
            = buffer.data() + (variant_idx * bytes_per_variant_);

        for (Eigen::Index sample_idx = 0; sample_idx < num_samples;
             ++sample_idx)
        {
            const Eigen::Index byte_index = sample_idx / 4;
            const Eigen::Index bit_offset = 2 * (sample_idx % 4);
            const uint8_t genotype_code
                = (variant_data[byte_index] >> bit_offset) & 0x03;

            raw_genotypes(sample_idx, variant_idx) = add_map_[genotype_code];
        }
    }

    // Reorder samples for each variant using the provided ID map
    Eigen::MatrixXd reordered_genotypes(id_map.size(), num_variants_in_chunk);

    for (Eigen::Index variant_idx = 0; variant_idx < num_variants_in_chunk;
         ++variant_idx)
    {
        reordered_genotypes.col(variant_idx)
            = reorder_genotypes(raw_genotypes.col(variant_idx), id_map);
    }

    return reordered_genotypes;
}

auto BedPipe::load() const -> std::expected<Eigen::MatrixXd, Error>
{
    const auto num_samples = static_cast<Eigen::Index>(load_sample_map_.size());

    if (num_samples == 0)
    {
        return std::unexpected(
            Error{ErrorCode::InvalidData, "No samples available for loading"});
    }

    // Use bulk reading for all variants
    return read_variants_bulk(0, num_variants(), load_sample_map_);
}

auto BedPipe::load(const std::unordered_map<std::string, Eigen::Index>& id_map)
    const -> std::expected<Eigen::MatrixXd, Error>
{
    if (id_map.empty())
    {
        return std::unexpected(
            Error{ErrorCode::InvalidData, "No samples available for loading"});
    }

    // Use bulk reading for all variants
    return read_variants_bulk(0, num_variants(), id_map);
}

auto BedPipe::load_chunk(Eigen::Index start_variant, Eigen::Index end_variant)
    const -> std::expected<Eigen::MatrixXd, Error>
{
    const auto num_samples = static_cast<Eigen::Index>(load_sample_map_.size());

    if (num_samples == 0)
    {
        return std::unexpected(
            Error{ErrorCode::InvalidData, "No samples available for loading"});
    }

    if (auto validation = validate_variant_index(start_variant); !validation)
    {
        return std::unexpected(validation.error());
    }
    if (auto validation = validate_variant_index(end_variant - 1); !validation)
    {
        return std::unexpected(validation.error());
    }

    if (start_variant >= end_variant)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidRange,
                "Invalid variant range for chunk loading"});
    }

    // Use bulk reading for the chunk
    return read_variants_bulk(start_variant, end_variant, load_sample_map_);
}

auto BedPipe::load_chunk(
    Eigen::Index start_variant,
    Eigen::Index end_variant,
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
    -> std::expected<Eigen::MatrixXd, Error>
{
    if (id_map.empty())
    {
        return std::unexpected(
            Error{ErrorCode::InvalidData, "No samples available for loading"});
    }

    if (auto validation = validate_variant_index(start_variant); !validation)
    {
        return std::unexpected(validation.error());
    }
    if (auto validation = validate_variant_index(end_variant - 1); !validation)
    {
        return std::unexpected(validation.error());
    }

    if (start_variant >= end_variant)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidRange,
                "Invalid variant range for chunk loading"});
    }

    // Use bulk reading for the chunk
    return read_variants_bulk(start_variant, end_variant, id_map);
}

auto BedPipe::intersect_samples(
    const std::unordered_set<std::string>& sample_ids)
    -> std::expected<void, Error>
{
    std::unordered_map<std::string, Eigen::Index> new_load_map;
    Eigen::Index new_index = 0;

    for (const auto& id : sample_ids)
    {
        if (raw_sample_map_.contains(id))
        {
            new_load_map[id] = new_index++;
        }
    }

    if (new_load_map.empty())
    {
        load_sample_map_.clear();

        auto logger = gelex::logging::get();
        logger->warn(
            "No matching samples found in intersection for BED file [{}]",
            bed_path_.string());

        return std::unexpected(
            Error{
                ErrorCode::InvalidData,
                "No matching samples found in intersection"});
    }

    auto logger = gelex::logging::get();
    logger->info(
        "Intersected {} samples from BED file [{}] with {} target samples, {} "
        "matches found",
        raw_sample_map_.size(),
        bed_path_.string(),
        sample_ids.size(),
        new_load_map.size());

    load_sample_map_ = std::move(new_load_map);
    return {};
}
}  // namespace gelex
