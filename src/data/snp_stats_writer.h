#pragma once

#include <expected>
#include <filesystem>
#include <fstream>
#include <vector>

#include <Eigen/Core>

#include "gelex/error.h"

namespace gelex::detail
{

class SnpStatsWriter
{
   public:
    static auto create(const std::filesystem::path& file_path)
        -> std::expected<SnpStatsWriter, Error>;

    SnpStatsWriter(const SnpStatsWriter&) = delete;
    SnpStatsWriter(SnpStatsWriter&&) noexcept = default;
    SnpStatsWriter& operator=(const SnpStatsWriter&) = delete;
    SnpStatsWriter& operator=(SnpStatsWriter&&) noexcept = default;
    ~SnpStatsWriter();

    auto write(
        Eigen::Index num_samples,
        Eigen::Index num_variants,
        size_t num_monomorphic,
        const std::vector<size_t>& monomorphic_indices,
        const std::vector<double>& means,
        const std::vector<double>& stddevs) -> std::expected<void, Error>;

   private:
    explicit SnpStatsWriter(
        std::ofstream&& file,
        std::filesystem::path&& file_path);
    std::ofstream file_;
    std::filesystem::path path_;
};

}  // namespace gelex::detail
