#pragma once

#include <array>
#include <expected>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader.h"
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"
#include "gelex/mio.h"

namespace gelex
{

struct VariantInstruction
{
    std::string id;
    Eigen::Index file_idx;
    bool reverse;
};

class BedPipe
{
   public:
    static auto create(
        const std::filesystem::path& bed_prefix,
        std::shared_ptr<SampleManager> sample_manager)
        -> std::expected<BedPipe, Error>;

    BedPipe(const BedPipe&) = delete;
    BedPipe& operator=(const BedPipe&) = delete;
    BedPipe(BedPipe&&) noexcept = default;
    BedPipe& operator=(BedPipe&&) noexcept = default;
    ~BedPipe() = default;

    auto load() const -> std::expected<Eigen::MatrixXd, Error>;
    auto load_chunk(Eigen::Index start_col, Eigen::Index end_col) const
        -> std::expected<Eigen::MatrixXd, Error>;

    void set_read_plan(std::vector<VariantInstruction> instructions);
    void reset_to_default();

    Eigen::Index num_samples() const;
    Eigen::Index num_variants() const;

    const std::vector<std::string>& snp_ids() const;

    static auto format_bed_path(std::string_view bed_path)
        -> std::expected<std::filesystem::path, Error>;

   private:
    BedPipe(
        mio::mmap_source&& mmap,
        std::unique_ptr<detail::BimLoader> bim_loader,
        std::shared_ptr<SampleManager> sample_manager,
        std::vector<Eigen::Index> sample_mapping,
        Eigen::Index raw_sample_count,
        Eigen::Index bytes_per_variant,
        std::filesystem::path bed_path);

    static auto validate_magic(const mio::mmap_source& mmap) -> bool;
    static auto calculate_bytes_per_variant(Eigen::Index num_samples)
        -> Eigen::Index;

    void decode_variant(
        const uint8_t* data_ptr,
        bool is_reverse,
        Eigen::Ref<Eigen::VectorXd> target_col) const;

    mio::mmap_source mmap_;
    std::unique_ptr<detail::BimLoader> bim_loader_;
    std::shared_ptr<SampleManager> sample_manager_;

    std::vector<Eigen::Index> raw_to_target_sample_idx_;

    std::vector<VariantInstruction> plan_;

    Eigen::Index raw_sample_count_;
    Eigen::Index bytes_per_variant_;
    std::filesystem::path bed_path_;
    bool use_custom_plan_;

    constexpr static std::array<double, 4> lut_vals_ = {2.0, 1.0, 1.0, 0.0};
    constexpr static std::array<double, 4> lut_vals_rev_ = {0.0, 1.0, 1.0, 2.0};
};

}  // namespace gelex
