#ifndef GELEX_DATA_GRM_LOADER_H
#define GELEX_DATA_GRM_LOADER_H

#include <cstddef>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include <mio.h>
#include <Eigen/Core>

#include "../src/types/freq_effect.h"

namespace gelex::detail
{

class GrmLoader
{
   public:
    // prefix: path without .grm.bin / .grm.id suffix
    explicit GrmLoader(const std::filesystem::path& prefix);

    GrmLoader(const GrmLoader&) = delete;
    GrmLoader(GrmLoader&&) noexcept = default;
    GrmLoader& operator=(const GrmLoader&) = delete;
    GrmLoader& operator=(GrmLoader&&) noexcept = default;
    ~GrmLoader() = default;

    // Load complete GRM matrix
    [[nodiscard]] auto load() const -> Eigen::MatrixXd;

    // Load GRM with filtering and reordering based on id_map
    // id_map: key = "FID_IID" format ID, value = target matrix row/col index
    // Throws InvalidInputException if any ID in id_map is not found in file
    [[nodiscard]] auto load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> Eigen::MatrixXd;

    [[nodiscard]] auto sample_ids() const noexcept
        -> const std::vector<std::string>&
    {
        return sample_ids_;
    }

    [[nodiscard]] auto num_samples() const noexcept -> Eigen::Index
    {
        return num_samples_;
    }

    [[nodiscard]] auto type() const noexcept -> freq::GrmType { return type_; }

   private:
    std::filesystem::path bin_path_;
    std::filesystem::path id_path_;
    mio::mmap_source mmap_;
    std::vector<std::string> sample_ids_;  // "FID_IID" format
    Eigen::Index num_samples_{};
    freq::GrmType type_;

    auto load_sample_ids() -> void;
    auto init_mmap() -> void;

    // Linear index in lower triangle for position (i, j) where i >= j
    [[nodiscard]] static auto lower_triangle_index(
        Eigen::Index i,
        Eigen::Index j) noexcept -> size_t
    {
        return (i * (i + 1) / 2) + j;
    }
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_GRM_LOADER_H
