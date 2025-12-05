#include "gelex/data/bed_pipe.h"

#include <omp.h>
#include <limits>
#include <system_error>

#include "../src/data/parser.h"
#include "gelex/exception.h"
#include "loader/fam_loader.h"

namespace gelex
{

consteval std::array<double, 4> generate_lut_entry(uint8_t byte, bool reverse)
{
    std::array<double, 4> entry{};

    // 00 -> Homozygote 1
    // 01 -> Missing
    // 10 -> Heterozygote
    // 11 -> Homozygote 2

    constexpr double nan = std::numeric_limits<double>::quiet_NaN();

    // IsReverse=false:
    // 00(0)->2.0, 01(1)->NaN, 10(2)->1.0, 11(3)->0.0
    constexpr std::array<double, 4> std_map = {2.0, nan, 1.0, 0.0};

    // IsReverse=true, swap A1 A2:
    // 00->0.0, 01->NaN, 10->1.0, 11->2.0
    constexpr std::array<double, 4> rev_map = {0.0, nan, 1.0, 2.0};

    const std::array<double, 4> map = reverse ? rev_map : std_map;

    for (int i = 0; i < 4; ++i)
    {
        entry[i] = map[(byte >> (2 * i)) & 0x03];
    }
    return entry;
}

template <bool Reverse>
consteval std::array<std::array<double, 4>, 256> generate_full_lut()
{
    std::array<std::array<double, 4>, 256> table{};
    for (size_t i = 0; i < 256; ++i)
    {
        table[i] = generate_lut_entry(static_cast<uint8_t>(i), Reverse);
    }
    return table;
}

const std::array<std::array<double, 4>, 256> BedPipe::kDecodeLut
    = generate_full_lut<false>();
const std::array<std::array<double, 4>, 256> BedPipe::kDecodeLutRev
    = generate_full_lut<true>();

BedPipe::BedPipe(
    const std::filesystem::path& bed_prefix,
    std::shared_ptr<SampleManager> sample_manager)
    : sample_manager_(std::move(sample_manager))
{
    if (!sample_manager_)
    {
        throw ArgumentValidationException("SampleManager cannot be null");
    }

    auto bed_path = bed_prefix;
    bed_path.replace_extension(".bed");
    auto bim_path = bed_prefix;
    bim_path.replace_extension(".bim");
    auto fam_path = bed_prefix;
    fam_path.replace_extension(".fam");

    bed_path_ = bed_path;

    num_raw_snps_
        = static_cast<Eigen::Index>(detail::count_total_lines(bim_path));

    std::unique_ptr<detail::FamLoader> fam_loader
        = std::make_unique<detail::FamLoader>(fam_path, false);

    const auto& raw_ids = fam_loader->ids();
    num_raw_samples_ = static_cast<Eigen::Index>(raw_ids.size());
    bytes_per_variant_ = (num_raw_samples_ + 3) / 4;

    const auto& target_map = sample_manager_->common_id_map();
    raw_to_target_sample_idx_.assign(num_raw_samples_, -1);

    Eigen::Index mapped_count = 0;
    bool sequential = true;

    for (Eigen::Index i = 0; i < num_raw_samples_; ++i)
    {
        if (auto it = target_map.find(raw_ids[i]); it != target_map.end())
        {
            raw_to_target_sample_idx_[i] = it->second;
            mapped_count++;

            if (it->second != i)
            {
                sequential = false;
            }
        }
        else
        {
            sequential = false;
        }
    }

    is_dense_mapping_ = sequential && (mapped_count == num_raw_samples_)
                        && (static_cast<size_t>(num_raw_samples_)
                            == sample_manager_->num_common_samples());

    std::error_code ec;
    mmap_.map(bed_path.string(), ec);
    if (ec)
    {
        throw FileOpenException(
            std::format("{}: failed to mmap bed file", bed_path.string()));
    }

    if (mmap_.size() < 3)
    {
        throw FileFormatException(
            std::format("{}: bed file too short", bed_path.string()));
    }
    if (mmap_[0] != 0x6C || mmap_[1] != 0x1B || mmap_[2] != 0x01)
    {
        throw FileFormatException(
            std::format("{}: invalid BED magic number", bed_path.string()));
    }

    size_t expected_size
        = 3 + (static_cast<size_t>(num_raw_snps_) * bytes_per_variant_);
    if (mmap_.size() < expected_size)
    {
        throw FileFormatException(
            std::format(
                "{}: bed file truncated. Expected {} bytes, got {}",
                bed_path.string(),
                expected_size,
                mmap_.size()));
    }
}
auto BedPipe::format_bed_path(std::string_view bed_path)
    -> std::filesystem::path
{
    std::filesystem::path bed
        = bed_path.contains(".bed") ? bed_path : std::string(bed_path) + ".bed";
    if (!std::filesystem::exists(bed))
    {
        throw FileNotFoundException(
            std::format("{}: file not found", bed.string()));
    }
    return bed;
}

void BedPipe::decode_variant_dense(
    const uint8_t* data_ptr,
    bool is_reverse,
    std::span<double> target_buf) const
{
    const auto& lut = is_reverse ? kDecodeLutRev : kDecodeLut;
    const Eigen::Index num_bytes = bytes_per_variant_;
    const Eigen::Index total_samples = num_raw_samples_;

    for (Eigen::Index i = 0; i < num_bytes; ++i)
    {
        const uint8_t byte = data_ptr[i];
        const auto& vals = lut[byte];  // 查表获取 4 个 double

        const Eigen::Index base_idx = i * 4;

        if (base_idx + 4 <= total_samples)
        {
            std::memcpy(&target_buf[base_idx], vals.data(), 4 * sizeof(double));
        }
        else
        {
            for (int k = 0; k < 4 && (base_idx + k) < total_samples; ++k)
            {
                target_buf[base_idx + k] = vals[k];
            }
        }
    }
}

void BedPipe::decode_variant_sparse(
    const uint8_t* data_ptr,
    bool is_reverse,
    std::span<double> target_buf) const
{
    const auto& lut = is_reverse ? kDecodeLutRev : kDecodeLut;
    const Eigen::Index num_bytes = bytes_per_variant_;

    for (Eigen::Index i = 0; i < num_bytes; ++i)
    {
        const uint8_t byte = data_ptr[i];
        const auto& vals = lut[byte];

        const Eigen::Index base_raw_idx = i * 4;

        if (Eigen::Index tidx = raw_to_target_sample_idx_[base_raw_idx];
            tidx != -1)
        {
            target_buf[tidx] = vals[0];
        }

        if (base_raw_idx + 1 < num_raw_samples_) [[likely]]
        {
            if (Eigen::Index tidx = raw_to_target_sample_idx_[base_raw_idx + 1];
                tidx != -1)
            {
                target_buf[tidx] = vals[1];
            }
        }

        if (base_raw_idx + 2 < num_raw_samples_) [[likely]]
        {
            if (Eigen::Index tidx = raw_to_target_sample_idx_[base_raw_idx + 2];
                tidx != -1)
            {
                target_buf[tidx] = vals[2];
            }
        }

        if (base_raw_idx + 3 < num_raw_samples_) [[likely]]
        {
            if (Eigen::Index tidx = raw_to_target_sample_idx_[base_raw_idx + 3];
                tidx != -1)
            {
                target_buf[tidx] = vals[3];
            }
        }
    }
}

Eigen::MatrixXd BedPipe::load_chunk(
    Eigen::Index start_col,
    Eigen::Index end_col) const
{
    const Eigen::Index max_cols = num_snps();

    if (start_col < 0 || end_col > max_cols || start_col >= end_col)
    {
        throw ColumnRangeException(
            std::format(
                "invalid chunk range: [{}, {}). Total SNPs: {}",
                start_col,
                end_col,
                max_cols));
    }

    const Eigen::Index num_output_rows = num_samples();
    const Eigen::Index num_output_cols = end_col - start_col;

    Eigen::MatrixXd result(num_output_rows, num_output_cols);

    if (!is_dense_mapping_)
    {
        result.setConstant(std::numeric_limits<double>::quiet_NaN());
    }

    const uint8_t* base_ptr
        = reinterpret_cast<const uint8_t*>(mmap_.data()) + 3;
    const size_t file_size = mmap_.size();

#pragma omp parallel for schedule(dynamic)
    for (Eigen::Index j = 0; j < num_output_cols; ++j)
    {
        const Eigen::Index current_col_idx = start_col + j;

        Eigen::Index target_matrix_col = j;
        size_t offset = 0;
        bool reverse = false;
        bool skip = false;

        if (use_custom_plan_)
        {
            const auto& instruction = plan_[current_col_idx];
            if (instruction.type == MatchType::skip)
            {
                skip = true;
            }
            else
            {
                target_matrix_col = instruction.target_col;
                offset = current_col_idx * bytes_per_variant_;
                reverse = (instruction.type == MatchType::reverse);
            }
        }
        else
        {
            offset = current_col_idx * bytes_per_variant_;
        }

        if (skip)
        {
            continue;
        }

        if (3 + offset + bytes_per_variant_ > file_size)
        {
            result.col(target_matrix_col)
                .setConstant(std::numeric_limits<double>::quiet_NaN());
            continue;
        }

        double* col_data_ptr = result.col(target_matrix_col).data();
        std::span<double> target_span(col_data_ptr, num_output_rows);

        const uint8_t* src_ptr = base_ptr + offset;

        if (is_dense_mapping_)
        {
            decode_variant_dense(src_ptr, reverse, target_span);
        }
        else
        {
            decode_variant_sparse(src_ptr, reverse, target_span);
        }
    }

    return result;
}

Eigen::MatrixXd BedPipe::load() const
{
    return load_chunk(0, num_snps());
}

void BedPipe::set_read_plan(MatchPlan&& instructions)
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

Eigen::Index BedPipe::num_snps() const
{
    if (use_custom_plan_)
    {
        return static_cast<Eigen::Index>(plan_.size());
    }
    return num_raw_snps_;
}

}  // namespace gelex
