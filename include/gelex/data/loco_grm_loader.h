#ifndef GELEX_DATA_LOCO_GRM_LOADER_H_
#define GELEX_DATA_LOCO_GRM_LOADER_H_

#include <filesystem>
#include <string>
#include <unordered_map>

#include <Eigen/Core>

namespace gelex
{
namespace detail
{
class GrmLoader;
}

class LocoGRMLoader
{
   public:
    LocoGRMLoader(
        const std::filesystem::path& whole_grm_prefix,
        const std::unordered_map<std::string, Eigen::Index>& id_map);

    LocoGRMLoader(const LocoGRMLoader&) = delete;
    LocoGRMLoader(LocoGRMLoader&&) noexcept = default;
    LocoGRMLoader& operator=(const LocoGRMLoader&) = delete;
    LocoGRMLoader& operator=(LocoGRMLoader&&) noexcept = default;
    ~LocoGRMLoader();

    /**
     * @brief Load the LOCO GRM for a specific chromosome.
     *
     * Formula: G_loco = (G_whole - G_i) / (K_whole - K_i)
     *
     * @param chr_grm_prefix Path prefix for the chromosome-specific GRM files.
     * @return Eigen::MatrixXd The calculated LOCO GRM.
     */
    [[nodiscard]] auto load_loco_grm(
        const std::filesystem::path& chr_grm_prefix,
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> Eigen::MatrixXd;

    [[nodiscard]] auto num_samples() const noexcept -> Eigen::Index;

   private:
    Eigen::MatrixXd
        g_whole_;       // Pre-filtered unnormalized whole matrix (Z * Z')
    double k_whole_{};  // K_whole
};

}  // namespace gelex

#endif  // GELEX_DATA_LOCO_GRM_LOADER_H_
