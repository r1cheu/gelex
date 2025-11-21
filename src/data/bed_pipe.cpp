#include "gelex/data/bed_pipe.h"

#include <omp.h>
#include <system_error>

namespace gelex
{

auto BedPipe::format_bed_path(std::string_view bed_path)
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

auto BedPipe::validate_magic(const mio::mmap_source& mmap) -> bool
{
    if (mmap.size() < 3)
    {
        return false;
    }
    return mmap[0] == 0x6C && mmap[1] == 0x1B && mmap[2] == 0x01;  // SNP-major
}

auto BedPipe::calculate_bytes_per_variant(Eigen::Index num_samples)
    -> Eigen::Index
{
    return (num_samples + 3) / 4;
}

auto BedPipe::create(
    const std::filesystem::path& bed_prefix,
    std::shared_ptr<SampleManager> sample_manager)
    -> std::expected<BedPipe, Error>
{
    auto bed_path = bed_prefix;
    bed_path.replace_extension(".bed");
    auto bim_path = bed_prefix;
    bim_path.replace_extension(".bim");
    auto fam_path = bed_prefix;
    fam_path.replace_extension(".fam");

    auto bim_loader = detail::BimLoader::create(bim_path);
    if (!bim_loader)
    {
        return std::unexpected(bim_loader.error());
    }

    auto fam_loader = detail::FamLoader::create(fam_path, false);
    if (!fam_loader)
    {
        return std::unexpected(fam_loader.error());
    }

    std::error_code err;
    mio::mmap_source mmap;
    mmap.map(bed_path.string(), err);
    if (err)
    {
        return std::unexpected(enrich_with_file_info(
            Error{
                ErrorCode::FileIOError,
                "Failed to memory-map BED file: " + err.message()},
            bed_path));
    }

    if (!validate_magic(mmap))
    {
        return std::unexpected(enrich_with_file_info(
            Error{ErrorCode::InvalidFile, "Invalid BED magic number or mode"},
            bed_path));
    }

    const auto& raw_ids = fam_loader->ids();
    const auto& target_map = sample_manager->common_id_map();

    const auto raw_count = static_cast<Eigen::Index>(raw_ids.size());
    std::vector<Eigen::Index> sample_mapping(raw_count, -1);

    for (Eigen::Index i = 0; i < raw_count; ++i)
    {
        const auto& id = raw_ids[i];
        if (auto it = target_map.find(id); it != target_map.end())
        {
            sample_mapping[i] = it->second;
        }
    }

    const Eigen::Index bytes_per_var = calculate_bytes_per_variant(raw_count);

    const auto num_snps = static_cast<Eigen::Index>(bim_loader->ids().size());
    const size_t expected_size = 3 + (num_snps * bytes_per_var);
    if (mmap.size() < expected_size)
    {
        return std::unexpected(enrich_with_file_info(
            Error{
                ErrorCode::InvalidFile,
                "BED file truncated or smaller than expected"},
            bed_path));
    }

    return BedPipe(
        std::move(mmap),
        std::make_unique<detail::BimLoader>(std::move(*bim_loader)),
        std::move(sample_manager),
        std::move(sample_mapping),
        raw_count,
        bytes_per_var,
        bed_path);
}

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------

BedPipe::BedPipe(
    mio::mmap_source&& mmap,
    std::unique_ptr<detail::BimLoader> bim_loader,
    std::shared_ptr<SampleManager> sample_manager,
    std::vector<Eigen::Index> sample_mapping,
    Eigen::Index raw_sample_count,
    Eigen::Index bytes_per_variant,
    std::filesystem::path bed_path)
    : mmap_(std::move(mmap)),
      bim_loader_(std::move(bim_loader)),
      sample_manager_(std::move(sample_manager)),
      raw_to_target_sample_idx_(std::move(sample_mapping)),
      raw_sample_count_(raw_sample_count),
      bytes_per_variant_(bytes_per_variant),
      bed_path_(std::move(bed_path)),
      use_custom_plan_(false)
{
}

void BedPipe::set_read_plan(std::vector<VariantInstruction> instructions)
{
    plan_ = std::move(instructions);
    use_custom_plan_ = true;
}

void BedPipe::reset_to_default()
{
    plan_.clear();
    use_custom_plan_ = false;
}

Eigen::Index BedPipe::num_samples() const
{
    return static_cast<Eigen::Index>(sample_manager_->num_common_samples());
}

Eigen::Index BedPipe::num_variants() const
{
    if (use_custom_plan_)
    {
        return static_cast<Eigen::Index>(plan_.size());
    }
    return static_cast<Eigen::Index>(bim_loader_->ids().size());
}

const std::vector<std::string>& BedPipe::snp_ids() const
{
    return bim_loader_->ids();
}

void BedPipe::decode_variant(
    const uint8_t* data_ptr,
    bool is_reverse,
    Eigen::Ref<Eigen::VectorXd> target_col) const
{
    const auto& lut = is_reverse ? lut_vals_rev_ : lut_vals_;

    const Eigen::Index raw_limit = raw_sample_count_;

    for (Eigen::Index i = 0; i < raw_limit; i += 4)
    {
        const uint8_t byte = data_ptr[i / 4];

        if (Eigen::Index target_idx = raw_to_target_sample_idx_[i];
            target_idx != -1)
        {
            target_col[target_idx] = lut[byte & 0x03];
        }

        if (i + 1 < raw_limit)
        {
            if (Eigen::Index target_idx = raw_to_target_sample_idx_[i + 1];
                target_idx != -1)
            {
                target_col[target_idx] = lut[(byte >> 2) & 0x03];
            }
        }

        if (i + 2 < raw_limit)
        {
            if (Eigen::Index target_idx = raw_to_target_sample_idx_[i + 2];
                target_idx != -1)
            {
                target_col[target_idx] = lut[(byte >> 4) & 0x03];
            }
        }

        if (i + 3 < raw_limit)
        {
            if (Eigen::Index target_idx = raw_to_target_sample_idx_[i + 3];
                target_idx != -1)
            {
                target_col[target_idx] = lut[(byte >> 6) & 0x03];
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Read Operations
// -----------------------------------------------------------------------------

auto BedPipe::load_chunk(Eigen::Index start_col, Eigen::Index end_col) const
    -> std::expected<Eigen::MatrixXd, Error>
{
    const Eigen::Index max_cols = num_variants();

    if (start_col < 0 || end_col > max_cols || start_col >= end_col)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidRange,
                std::format(
                    "Invalid chunk range: [{}, {}). Total: {}",
                    start_col,
                    end_col,
                    max_cols)});
    }

    const Eigen::Index num_output_rows = num_samples();
    const Eigen::Index num_output_cols = end_col - start_col;

    Eigen::MatrixXd result(num_output_rows, num_output_cols);

    const uint8_t* base_ptr
        = reinterpret_cast<const uint8_t*>(mmap_.data()) + 3;
    const size_t file_size = mmap_.size();

#pragma omp parallel for schedule(dynamic)
    for (Eigen::Index j = 0; j < num_output_cols; ++j)
    {
        const Eigen::Index current_col_idx = start_col + j;

        size_t offset = 0;
        bool reverse = false;

        if (use_custom_plan_)
        {
            const auto& instruction = plan_[current_col_idx];
            offset = instruction.file_idx * bytes_per_variant_;
            reverse = instruction.reverse;
        }
        else
        {
            offset = current_col_idx * bytes_per_variant_;
            reverse = false;
        }

        if (3 + offset + bytes_per_variant_ > file_size)
        {
            result.col(j).setConstant(-9.0);
            continue;
        }

        decode_variant(base_ptr + offset, reverse, result.col(j));
    }

    return result;
}

auto BedPipe::load() const -> std::expected<Eigen::MatrixXd, Error>
{
    return load_chunk(0, num_variants());
}

}  // namespace gelex
