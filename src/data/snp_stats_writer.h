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
    explicit SnpStatsWriter(std::filesystem::path file_path);

    SnpStatsWriter(const SnpStatsWriter&) = delete;
    SnpStatsWriter(SnpStatsWriter&&) noexcept = default;
    SnpStatsWriter& operator=(const SnpStatsWriter&) = delete;
    SnpStatsWriter& operator=(SnpStatsWriter&&) noexcept = default;
    ~SnpStatsWriter();

    auto open() -> std::expected<void, Error>;
    auto close() -> std::expected<void, Error>;

    auto write_all(
        Eigen::Index num_samples,
        Eigen::Index num_variants,
        size_t num_monomorphic,
        const std::vector<size_t>& monomorphic_indices,
        const std::vector<double>& means,
        const std::vector<double>& stddevs) -> std::expected<void, Error>;

    const std::filesystem::path& file_path() const noexcept
    {
        return file_path_;
    }
    bool is_open() const noexcept { return file_.is_open(); }

   private:
    std::filesystem::path file_path_;
    std::ofstream file_;
};

}  // namespace gelex::detail
