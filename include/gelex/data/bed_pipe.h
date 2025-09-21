#pragma once

#include <expected>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>

#include "../src/data/loader.h"
#include "Eigen/Core"
#include "gelex/error.h"

namespace gelex
{

class BedPipe
{
   public:
    static auto create(
        const std::filesystem::path& prefix,
        bool iid_only = false) -> std::expected<BedPipe, Error>;

    BedPipe(const BedPipe&) = delete;
    BedPipe(BedPipe&&) noexcept = default;
    BedPipe& operator=(const BedPipe&) = delete;
    BedPipe& operator=(BedPipe&&) noexcept = default;
    ~BedPipe() = default;

    auto get_genotypes(Eigen::Index variant_index) const
        -> std::expected<Eigen::VectorXd, Error>;
    auto get_sample_genotypes(Eigen::Index sample_index) const
        -> std::expected<Eigen::VectorXd, Error>;
    auto get_genotype(Eigen::Index variant_index, Eigen::Index sample_index)
        const -> std::expected<double, Error>;

    Eigen::Index num_variants() const
    {
        return static_cast<Eigen::Index>(bim_loader_->snp_ids().size());
    }
    Eigen::Index sample_size() const
    {
        return static_cast<Eigen::Index>(fam_loader_->sample_map().size());
    }

    const std::vector<std::string>& snp_ids() const
    {
        return bim_loader_->snp_ids();
    }
    const std::unordered_set<std::string>& sample_ids() const
    {
        return fam_loader_->sample_ids();
    }
    const std::unordered_map<std::string, Eigen::Index>& sample_map() const
    {
        return fam_loader_->sample_map();
    }

    auto load(
        const std::optional<std::unordered_map<std::string, Eigen::Index>>&
            id_map
        = std::nullopt) const -> std::expected<Eigen::MatrixXd, Error>;
    auto load_chunk(
        Eigen::Index start_variant,
        Eigen::Index end_variant,
        const std::optional<std::unordered_map<std::string, Eigen::Index>>&
            id_map
        = std::nullopt) const -> std::expected<Eigen::MatrixXd, Error>;

   private:
    BedPipe(
        std::ifstream&& file_stream,
        std::unique_ptr<detail::FamLoader> fam_loader,
        std::unique_ptr<detail::BimLoader> bim_loader,
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

    auto validate_id_map(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> std::expected<void, Error>;

    Eigen::VectorXd reorder_genotypes(
        const Eigen::VectorXd& raw_genotypes,
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    auto read_variants_bulk(
        Eigen::Index start_variant,
        Eigen::Index end_variant,
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> std::expected<Eigen::MatrixXd, Error>;

    mutable std::ifstream file_stream_;
    std::unique_ptr<detail::FamLoader> fam_loader_;
    std::unique_ptr<detail::BimLoader> bim_loader_;
    Eigen::Index bytes_per_variant_;
    std::filesystem::path bed_path_;

    constexpr static std::array<double, 4> encode_map_{0, 1, 1, 2};
};

}  // namespace gelex
