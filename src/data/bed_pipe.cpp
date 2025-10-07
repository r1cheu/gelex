#include "gelex/data/bed_pipe.h"

#include <expected>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader.h"
#include "gelex/error.h"
#include "gelex/logger.h"

namespace gelex
{

auto valid_bed(std::string_view bed_path)
    -> std::expected<std::filesystem::path, Error>
{
    std::filesystem::path bed
        = bed_path.contains(".bed") ? bed_path : std::string(bed_path) + ".bed";
    if (!std::filesystem::exists(bed))
    {
        return std::unexpected(
            Error{
                .code = ErrorCode::FileNotFound,
                .message = std::format("bed file [{}] not found", bed.string())}

        );
    }
    return bed;
}

BedPipe::BedPipe(
    std::ifstream&& file_stream,
    std::unique_ptr<detail::BimLoader> bim_loader,
    std::shared_ptr<SampleManager> sample_manager,
    const Eigen::Index bytes_per_variant,
    std::filesystem::path bed_path)
    : file_stream_(std::move(file_stream)),
      bim_loader_(std::move(bim_loader)),
      sample_manager_(std::move(sample_manager)),
      bytes_per_variant_(bytes_per_variant),
      bed_path_(std::move(bed_path))
{
    // File is already validated and opened by create()
    file_stream_.seekg(3);  // Skip the 3-byte magic number
}

auto BedPipe::create(
    const std::filesystem::path& bed_path,
    std::shared_ptr<SampleManager> sample_manager)
    -> std::expected<BedPipe, Error>
{
    // Create bim file path
    std::filesystem::path bim_path = bed_path;
    bim_path.replace_extension(".bim");

    // Load bim file
    auto bim_loader = detail::BimLoader::create(bim_path);
    if (!bim_loader)
    {
        return std::unexpected(bim_loader.error());
    }

    auto file = detail::open_file<std::ifstream>(
        bed_path, std::ios_base::in | std::ios_base::binary);
    if (!file)
    {
        return std::unexpected(file.error());
    }

    if (auto validation = validate_bed_file(*file, bed_path); !validation)
    {
        return std::unexpected(validation.error());
    }

    const auto num_samples
        = static_cast<Eigen::Index>(sample_manager->num_genotyped_samples());
    const Eigen::Index bytes_per_variant
        = calculate_bytes_per_variant(num_samples);

    auto logger = gelex::logging::get();
    logger->info(
        "Loaded BED file with {} samples and {} variants",
        sample_manager->num_common_samples(),
        bim_loader->ids().size());

    return BedPipe(
        std::move(*file),
        std::make_unique<detail::BimLoader>(std::move(*bim_loader)),
        std::move(sample_manager),
        bytes_per_variant,
        bed_path);
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
            path));
    }

    if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01)
    {
        return std::unexpected(enrich_with_file_info(
            Error{ErrorCode::InvalidFile, "Invalid BED file magic number"},
            path));
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
    const Eigen::VectorXd& raw_genotypes) const
{
    const auto& id_map = sample_manager_->common_id_map();
    const auto& genotyped_sample_map = sample_manager_->genotyped_sample_map();

    Eigen::VectorXd reordered(id_map.size());

    for (const auto& [sample_id, load_index] : id_map)
    {
        const Eigen::Index raw_index = genotyped_sample_map.at(sample_id);
        reordered(load_index) = raw_genotypes(raw_index);
    }

    return reordered;
}

auto BedPipe::get_genotypes(Eigen::Index variant_index) const
    -> std::expected<Eigen::VectorXd, Error>
{
    return get_genotypes_impl(variant_index, false);
}

auto BedPipe::get_genotypes_impl(Eigen::Index variant_index, bool dominant)
    const -> std::expected<Eigen::VectorXd, Error>
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
            bed_path_));
    }

    const auto num_samples
        = static_cast<Eigen::Index>(sample_manager_->num_genotyped_samples());
    Eigen::VectorXd genotypes(num_samples);

    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        const Eigen::Index byte_index = i / 4;
        const Eigen::Index bit_offset = 2 * (i % 4);
        const uint8_t genotype_code = (buffer[byte_index] >> bit_offset) & 0x03;
        genotypes(i)
            = dominant ? dom_map_[genotype_code] : add_map_[genotype_code];
    }

    return reorder_genotypes(genotypes);
}

auto BedPipe::get_genotype_impl(
    Eigen::Index variant_index,
    Eigen::Index sample_index,
    bool dominant) const -> std::expected<double, Error>
{
    const auto num_samples
        = static_cast<Eigen::Index>(sample_manager_->num_genotyped_samples());

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
            bed_path_));
    }

    const Eigen::Index bit_offset = 2 * (sample_index % 4);
    const uint8_t genotype_code = (byte >> bit_offset) & 0x03;

    return dominant ? dom_map_[genotype_code] : add_map_[genotype_code];
}

auto BedPipe::get_genotype(
    Eigen::Index variant_index,
    Eigen::Index sample_index) const -> std::expected<double, Error>
{
    return get_genotype_impl(variant_index, sample_index, false);
}

auto BedPipe::get_sample_genotypes(Eigen::Index sample_index) const
    -> std::expected<Eigen::VectorXd, Error>
{
    if (sample_index >= sample_size())
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

auto BedPipe::get_dominant_genotypes(Eigen::Index variant_index) const
    -> std::expected<Eigen::VectorXd, Error>
{
    return get_genotypes_impl(variant_index, true);
}

auto BedPipe::get_dominant_genotype(
    Eigen::Index variant_index,
    Eigen::Index sample_index) const -> std::expected<double, Error>
{
    return get_genotype_impl(variant_index, sample_index, true);
}

auto BedPipe::get_sample_dominant_genotypes(Eigen::Index sample_index) const
    -> std::expected<Eigen::VectorXd, Error>
{
    if (sample_index >= sample_size())
    {
        return std::unexpected(
            Error{ErrorCode::InvalidRange, "Sample index out of range"});
    }

    Eigen::VectorXd genotypes(num_variants());

    for (Eigen::Index i = 0; i < num_variants(); ++i)
    {
        auto result = get_dominant_genotype(i, sample_index);
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
    bool dominant) const -> std::expected<Eigen::MatrixXd, Error>
{
    const auto num_samples
        = static_cast<Eigen::Index>(sample_manager_->num_genotyped_samples());
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
            bed_path_));
    }

    // Parse all variants from the buffer
    Eigen::MatrixXd genotypes(num_samples, num_variants_in_chunk);

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

            genotypes(sample_idx, variant_idx)
                = dominant ? dom_map_[genotype_code] : add_map_[genotype_code];
        }
    }

    // Always reorder genotypes to match the intersected sample order
    Eigen::MatrixXd reordered_genotypes(
        sample_manager_->num_common_samples(), num_variants_in_chunk);

    for (Eigen::Index variant_idx = 0; variant_idx < num_variants_in_chunk;
         ++variant_idx)
    {
        reordered_genotypes.col(variant_idx)
            = reorder_genotypes(genotypes.col(variant_idx));
    }

    return reordered_genotypes;
}

auto BedPipe::load(bool dominant) const -> std::expected<Eigen::MatrixXd, Error>
{
    // Use bulk reading for all variants
    return read_variants_bulk(0, num_variants(), dominant);
}

auto BedPipe::load_chunk(
    Eigen::Index start_variant,
    Eigen::Index end_variant,
    bool dominant) const -> std::expected<Eigen::MatrixXd, Error>
{
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
    return read_variants_bulk(start_variant, end_variant, dominant);
}

}  // namespace gelex
