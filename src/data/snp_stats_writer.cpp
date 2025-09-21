#include "snp_stats_writer.h"

#include <expected>
#include <filesystem>
#include <fstream>
#include <vector>

#include <Eigen/Core>

#include "data/loader.h"
#include "gelex/error.h"

namespace gelex::detail
{

using Eigen::Index;

SnpStatsWriter::SnpStatsWriter(
    std::ofstream&& file,
    std::filesystem::path&& file_path)
    : file_(std::move(file)), path_(std::move(file_path))
{
}

auto SnpStatsWriter::create(const std::filesystem::path& file_path)
    -> std::expected<SnpStatsWriter, Error>
{
    auto file = open_file<std::ofstream>(file_path, std::ios::binary);

    if (!file)
    {
        return std::unexpected(file.error());
    }
    return SnpStatsWriter(std::move(*file), std::move(file_path.string()));
}

auto SnpStatsWriter::write(
    Index num_samples,
    Index num_variants,
    size_t num_monomorphic,
    const std::vector<size_t>& monomorphic_indices,
    const std::vector<double>& means,
    const std::vector<double>& stddevs) -> std::expected<void, Error>
{
    if (!file_.is_open())
    {
        return std::unexpected(
            Error{ErrorCode::FileIOError, "Stats file is not open"});
    }

    file_.write(
        reinterpret_cast<const char*>(&num_samples), sizeof(num_samples));
    file_.write(
        reinterpret_cast<const char*>(&num_variants), sizeof(num_variants));
    file_.write(
        reinterpret_cast<const char*>(&num_monomorphic),
        sizeof(num_monomorphic));

    if (!file_.good())
    {
        return std::unexpected(enrich_with_file_info(
            Error{
                ErrorCode::FileIOError, "Failed to write header to stats file"},
            path_));
    }

    if (!monomorphic_indices.empty())
    {
        file_.write(
            reinterpret_cast<const char*>(monomorphic_indices.data()),
            monomorphic_indices.size() * sizeof(size_t));
    }

    if (!means.empty())
    {
        file_.write(
            reinterpret_cast<const char*>(means.data()),
            means.size() * sizeof(double));
    }

    if (!stddevs.empty())
    {
        file_.write(
            reinterpret_cast<const char*>(stddevs.data()),
            stddevs.size() * sizeof(double));
    }

    if (!file_.good())
    {
        return std::unexpected(enrich_with_file_info(
            Error{ErrorCode::FileIOError, "Failed to write data to stats file"},
            path_));
    }

    return {};
}

}  // namespace gelex::detail
