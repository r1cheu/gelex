#include "snp_stats_writer.h"

#include <expected>
#include <filesystem>
#include <fstream>
#include <vector>

#include <Eigen/Core>

#include "gelex/error.h"

namespace gelex::detail
{

using Eigen::Index;

SnpStatsWriter::SnpStatsWriter(std::filesystem::path file_path)
    : file_path_(std::move(file_path))
{
}

SnpStatsWriter::~SnpStatsWriter()
{
    if (file_.is_open())
    {
        file_.close();
    }
}

auto SnpStatsWriter::open() -> std::expected<void, Error>
{
    file_.open(file_path_, std::ios::binary);
    if (!file_.is_open())
    {
        return std::unexpected(enrich_with_file_info(
            Error{ErrorCode::FileIOError, "Failed to open stats file"},
            file_path_.string()));
    }

    return {};
}

auto SnpStatsWriter::close() -> std::expected<void, Error>
{
    if (file_.is_open())
    {
        file_.close();
    }

    if (file_.fail())
    {
        return std::unexpected(enrich_with_file_info(
            Error{ErrorCode::FileIOError, "Failed to close stats file"},
            file_path_.string()));
    }

    return {};
}

auto SnpStatsWriter::write_all(
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

    // Write header information
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
            file_path_.string()));
    }

    // Write monomorphic indices
    if (!monomorphic_indices.empty())
    {
        file_.write(
            reinterpret_cast<const char*>(monomorphic_indices.data()),
            monomorphic_indices.size() * sizeof(size_t));
    }

    // Write means
    if (!means.empty())
    {
        file_.write(
            reinterpret_cast<const char*>(means.data()),
            means.size() * sizeof(double));
    }

    // Write standard deviations
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
            file_path_.string()));
    }

    return {};
}

}  // namespace gelex::detail
