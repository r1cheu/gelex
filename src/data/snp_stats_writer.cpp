#include "snp_stats_writer.h"

#include <cstdint>
#include <expected>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
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
    return SnpStatsWriter(std::move(*file), file_path.string());
}

auto SnpStatsWriter::write(
    int64_t num_samples,
    int64_t num_variants,
    int64_t num_monomorphic,
    const std::vector<int64_t>& monomorphic_indices,
    const std::vector<double>& means,
    const std::vector<double>& stddevs) -> std::expected<void, Error>
{
    if (!file_.is_open())
    {
        return std::unexpected(
            Error{ErrorCode::FileIOError, "Stats file is not open"});
    }

    if (num_monomorphic != static_cast<int64_t>(monomorphic_indices.size()))
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidArgument,
                "num_monomorphic must equal monomorphic_indices.size()"});
    }

    if (means.size() != static_cast<size_t>(num_variants)
        || stddevs.size() != static_cast<size_t>(num_variants))
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidArgument,
                "means.size() and stddevs.size() must equal num_variants"});
    }

    if (means.empty() || stddevs.empty())
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidArgument,
                "means and stddevs vectors cannot be empty"});
    }
    constexpr size_t int64_t_size = sizeof(int64_t);

    file_.write(reinterpret_cast<const char*>(&num_samples), int64_t_size);
    file_.write(reinterpret_cast<const char*>(&num_variants), int64_t_size);
    file_.write(reinterpret_cast<const char*>(&num_monomorphic), int64_t_size);

    if (!file_.good())
    {
        return std::unexpected(enrich_with_file_info(
            Error{
                ErrorCode::FileIOError, "Failed to write header to stats file"},
            path_));
    }

    file_.write(
        reinterpret_cast<const char*>(monomorphic_indices.data()),
        static_cast<std::streamsize>(
            monomorphic_indices.size() * int64_t_size));

    file_.write(
        reinterpret_cast<const char*>(means.data()),
        static_cast<std::streamsize>(means.size() * sizeof(double)));

    file_.write(
        reinterpret_cast<const char*>(stddevs.data()),
        static_cast<std::streamsize>(stddevs.size() * sizeof(double)));

    if (!file_.good())
    {
        return std::unexpected(enrich_with_file_info(
            Error{ErrorCode::FileIOError, "Failed to write data to stats file"},
            path_));
    }

    // Flush to ensure data is written to disk
    file_.flush();

    if (!file_.good())
    {
        return std::unexpected(enrich_with_file_info(
            Error{ErrorCode::FileIOError, "Failed to flush stats file"},
            path_));
    }

    return {};
}

}  // namespace gelex::detail
