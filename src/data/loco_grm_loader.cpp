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
    // Load and filter the whole GRM once during construction.
    // load_unnormalized(id_map) returns (X_w * X_w') filtered and reordered.
    g_whole_ = whole_loader.load_unnormalized(id_map);
    // Compute trace after loading and save for LOCO calculation
    trace_whole_ = g_whole_.trace();
    k_whole_ = trace_whole_ / static_cast<double>(g_whole_.rows());
}

LocoGRMLoader::~LocoGRMLoader() = default;

auto LocoGRMLoader::load_loco_grm(
    const std::filesystem::path& chr_grm_prefix,
    const std::unordered_map<std::string, Eigen::Index>& id_map,
    Eigen::MatrixXd& target) const -> void
{
    std::filesystem::path bin_path = chr_grm_prefix.string() + ".bin";
    if (!std::filesystem::exists(bin_path))
    {
        throw InvalidInputException(
            std::format(
                "LOCO error: GRM file not found: {}", bin_path.string()));
    }

    detail::GrmLoader chr_loader(chr_grm_prefix);

    // Load chromosome GRM filtered by the SAME id_map to ensure alignment.
    // Use the mutable buffer to avoid reallocations.
    chr_loader.load_unnormalized(id_map, g_chr_buffer_);

    // Compute k_i from the loaded chromosome GRM trace
    double trace_i = g_chr_buffer_.trace();
    double k_i = trace_i / static_cast<double>(g_chr_buffer_.rows());
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

    // Both g_whole_ and g_chr_buffer_ are now unnormalized (X*X') and filtered.
    target = (g_whole_ - g_chr_buffer_) / k_loco;
}

auto LocoGRMLoader::load_loco_grm(
    const std::filesystem::path& chr_grm_prefix,
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
    -> Eigen::MatrixXd
{
    Eigen::MatrixXd target;
    load_loco_grm(chr_grm_prefix, id_map, target);
    return target;
}

auto LocoGRMLoader::num_samples() const noexcept -> Eigen::Index
{
    return g_whole_.rows();
}

}  // namespace gelex
