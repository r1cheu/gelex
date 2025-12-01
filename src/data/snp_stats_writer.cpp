#include "snp_stats_writer.h"

#include <cassert>
#include <format>

#include "gelex/exception.h"
#include "parser.h"

namespace gelex::detail
{

SnpStatsWriter::SnpStatsWriter(const std::filesystem::path& file_path)
    : path_(file_path), io_buffer_(kDefaultBufferSize)
{
    file_ = detail::open_file<std::ofstream>(
        path_, std::ios::binary | std::ios::trunc, io_buffer_);
}

void SnpStatsWriter::write(
    int64_t num_samples,
    std::span<const int64_t> monomorphic_indices,
    std::span<const double> means,
    std::span<const double> stddevs)
{
    if (!file_.is_open())
    {
        throw FileOpenException(
            enrich_with_file_info("Stats file is not open", path_));
    }

    if (means.size() != stddevs.size())
    {
        throw ArgumentValidationException(
            std::format(
                "Size mismatch: means ({}) and stddevs ({}) must have the same "
                "length.",
                means.size(),
                stddevs.size()));
    }

    if (means.empty())
    {
        throw ArgumentValidationException("means and stddevs cannot be empty");
    }

    const auto num_variants = static_cast<int64_t>(means.size());
    const auto num_monomorphic
        = static_cast<int64_t>(monomorphic_indices.size());

    const int64_t header[] = {num_samples, num_variants, num_monomorphic};
    file_.write(reinterpret_cast<const char*>(header), sizeof(header));

    if (!file_.good())
    {
        throw FileOpenException(enrich_with_file_info(
            "Failed to write header to stats file", path_));
    }

    if (!monomorphic_indices.empty())
    {
        file_.write(
            reinterpret_cast<const char*>(monomorphic_indices.data()),
            static_cast<std::streamsize>(monomorphic_indices.size_bytes()));
    }

    file_.write(
        reinterpret_cast<const char*>(means.data()),
        static_cast<std::streamsize>(means.size_bytes()));

    file_.write(
        reinterpret_cast<const char*>(stddevs.data()),
        static_cast<std::streamsize>(stddevs.size_bytes()));

    if (!file_.good())
    {
        throw FileOpenException(
            enrich_with_file_info("Failed to write data to stats file", path_));
    }

    file_.flush();

    if (!file_.good())
    {
        throw FileOpenException(
            enrich_with_file_info("Failed to flush stats file", path_));
    }
}

}  // namespace gelex::detail
