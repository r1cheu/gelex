#include "gelex/data/loco_grm_loader.h"

#include <format>

#include "gelex/exception.h"
#include "grm_loader.h"

namespace gelex
{

LocoGRMLoader::LocoGRMLoader(
    const std::filesystem::path& whole_grm_prefix,
    const std::unordered_map<std::string, Eigen::Index>& id_map)
{
    detail::GrmLoader whole_loader(whole_grm_prefix);
    k_whole_ = whole_loader.denominator();
    // Load and filter the whole GRM once during construction.
    // load_unnormalized(id_map) returns (X_w * X_w') filtered and reordered.
    g_whole_ = whole_loader.load_unnormalized(id_map);
}

LocoGRMLoader::~LocoGRMLoader() = default;

auto LocoGRMLoader::load_loco_grm(
    const std::filesystem::path& chr_grm_prefix,
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
    -> Eigen::MatrixXd
{
    detail::GrmLoader chr_loader(chr_grm_prefix);

    // Load chromosome GRM filtered by the SAME id_map to ensure alignment.
    Eigen::MatrixXd g_i = chr_loader.load_unnormalized(id_map);

    double k_i = chr_loader.denominator();
    double k_loco = k_whole_ - k_i;
    if (k_loco <= 0)
    {
        throw InvalidInputException(
            std::format(
                "LOCO error: Chromosome GRM denominator ({}) is greater than "
                "or equal to Whole GRM denominator ({})",
                k_i,
                k_whole_));
    }

    // Both g_whole_ and g_i are now unnormalized (X*X') and filtered.
    return (g_whole_ - g_i) / k_loco;
}

auto LocoGRMLoader::num_samples() const noexcept -> Eigen::Index
{
    return g_whole_.rows();
}

}  // namespace gelex
