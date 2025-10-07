#pragma once

#include <expected>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "../src/data/loader.h"
#include "Eigen/Core"
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"

namespace gelex
{

auto valid_bed(std::string_view bed_path)
    -> std::expected<std::filesystem::path, Error>;

class BedPipe
{
   public:
    static auto create(
        const std::filesystem::path& bed_path,
        std::shared_ptr<SampleManager> sample_manager)
        -> std::expected<BedPipe, Error>;

    BedPipe(const BedPipe&) = delete;
    BedPipe(BedPipe&&) noexcept = default;
    BedPipe& operator=(const BedPipe&) = delete;
    BedPipe& operator=(BedPipe&&) noexcept = default;
    ~BedPipe() = default;

    auto get_genotypes(Eigen::Index variant_index) const
        -> std::expected<Eigen::VectorXd, Error>;
    auto get_dominant_genotypes(Eigen::Index variant_index) const
        -> std::expected<Eigen::VectorXd, Error>;
    auto get_sample_genotypes(Eigen::Index sample_index) const
        -> std::expected<Eigen::VectorXd, Error>;
    auto get_sample_dominant_genotypes(Eigen::Index sample_index) const
        -> std::expected<Eigen::VectorXd, Error>;
    auto get_genotype(Eigen::Index variant_index, Eigen::Index sample_index)
        const -> std::expected<double, Error>;
    auto get_dominant_genotype(
        Eigen::Index variant_index,
        Eigen::Index sample_index) const -> std::expected<double, Error>;

    Eigen::Index num_variants() const
    {
        return static_cast<Eigen::Index>(bim_loader_->ids().size());
    }
    Eigen::Index sample_size() const
    {
        return static_cast<Eigen::Index>(sample_manager_->num_common_samples());
    }

    const std::vector<std::string>& snp_ids() const
    {
        return bim_loader_->ids();
    }

    auto load(bool dominant = false) const
        -> std::expected<Eigen::MatrixXd, Error>;
    auto load_chunk(
        Eigen::Index start_variant,
        Eigen::Index end_variant,
        bool dominant = false) const -> std::expected<Eigen::MatrixXd, Error>;

   private:
    BedPipe(
        std::ifstream&& file_stream,
        std::unique_ptr<detail::BimLoader> bim_loader,
        std::shared_ptr<SampleManager> sample_manager,
        Eigen::Index bytes_per_variant,
        std::filesystem::path bed_path);

    static auto validate_bed_file(
        std::ifstream& file,
        const std::filesystem::path& path) -> std::expected<void, Error>;
    static auto calculate_bytes_per_variant(Eigen::Index num_samples)
        -> Eigen::Index;
    Eigen::Index calculate_offset(Eigen::Index variant_index) const;
    auto validate_variant_index(Eigen::Index variant_index) const
        -> std::expected<void, Error>;

    Eigen::VectorXd reorder_genotypes(
        const Eigen::VectorXd& raw_genotypes) const;

    auto read_variants_bulk(
        Eigen::Index start_variant,
        Eigen::Index end_variant,
        bool dominant = false) const -> std::expected<Eigen::MatrixXd, Error>;

    auto get_genotypes_impl(Eigen::Index variant_index, bool dominant) const
        -> std::expected<Eigen::VectorXd, Error>;
    auto get_genotype_impl(
        Eigen::Index variant_index,
        Eigen::Index sample_index,
        bool dominant) const -> std::expected<double, Error>;

    mutable std::ifstream file_stream_;
    std::unique_ptr<detail::BimLoader> bim_loader_;
    std::shared_ptr<SampleManager> sample_manager_;
    Eigen::Index bytes_per_variant_;
    std::filesystem::path bed_path_;

    constexpr static std::array<double, 4> add_map_{2, 1, 1, 0};
    constexpr static std::array<double, 4> dom_map_{0, 1, 1, 0};
};

}  // namespace gelex
