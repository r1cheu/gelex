#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <span>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/sample_manager.h"
#include "mio.h"

#include "../src/predictor/snp_matcher.h"

namespace gelex
{

class BedPipe
{
   public:
    BedPipe(
        const std::filesystem::path& bed_prefix,
        std::shared_ptr<SampleManager> sample_manager);

    BedPipe(const BedPipe&) = delete;
    BedPipe& operator=(const BedPipe&) = delete;
    BedPipe(BedPipe&&) noexcept = default;
    BedPipe& operator=(BedPipe&&) noexcept = default;
    ~BedPipe() = default;

    [[nodiscard]] Eigen::MatrixXd load() const;

    [[nodiscard]] Eigen::MatrixXd load_chunk(
        Eigen::Index start_col,
        Eigen::Index end_col) const;

    void set_read_plan(MatchPlan&& instructions);
    void reset_to_default();

    [[nodiscard]] Eigen::Index num_samples() const;
    [[nodiscard]] Eigen::Index num_snps() const;

    static auto format_bed_path(std::string_view bed_path)
        -> std::filesystem::path;

   private:
    void decode_variant_dense(
        const uint8_t* data_ptr,
        bool is_reverse,
        std::span<double> target_buf) const;

    void decode_variant_sparse(
        const uint8_t* data_ptr,
        bool is_reverse,
        std::span<double> target_buf) const;

    mio::mmap_source mmap_;
    std::shared_ptr<SampleManager> sample_manager_;

    std::vector<Eigen::Index> raw_to_target_sample_idx_;

    bool is_dense_mapping_ = false;

    MatchPlan plan_;
    bool use_custom_plan_ = false;

    Eigen::Index num_raw_samples_ = 0;
    Eigen::Index num_raw_snps_ = 0;
    Eigen::Index bytes_per_variant_ = 0;
    std::filesystem::path bed_path_;

    static constexpr size_t kLutSize = 256;
    using LutEntry = std::array<double, 4>;

    static const std::array<LutEntry, kLutSize> kDecodeLut;
    static const std::array<LutEntry, kLutSize> kDecodeLutRev;
};

}  // namespace gelex
