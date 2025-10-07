#include "gelex/data/genotype_mmap.h"

#include <cstdint>
#include <filesystem>
#include <fstream>

#include <fmt/ranges.h>
#include <gelex/mio.h>
#include <Eigen/Core>

#include "loader.h"

namespace gelex
{
using std::ifstream;

auto GenotypeMap::create(const std::filesystem::path& bin_file)
    -> std::expected<GenotypeMap, Error>
{
    auto snp_stats = bin_file;
    snp_stats.replace_extension(".snpstats");

    auto meta_stream_result = detail::open_file<ifstream>(
        snp_stats, std::ios::in | std::ios::binary);

    if (!meta_stream_result)
    {
        return std::unexpected(meta_stream_result.error());
    }

    auto& meta_stream = *meta_stream_result;

    int64_t rows = 0;
    int64_t cols = 0;
    int64_t num_mono_snp = 0;
    meta_stream.read(reinterpret_cast<char*>(&rows), sizeof(int64_t));
    if (!meta_stream)
    {
        return std::unexpected(
            Error{
                ErrorCode::FileIOError,
                std::format(
                    "Failed to read number of rows from metadata file: {}",
                    snp_stats.string())});
    }
    meta_stream.read(reinterpret_cast<char*>(&cols), sizeof(int64_t));
    if (!meta_stream)
    {
        return std::unexpected(
            Error{
                ErrorCode::FileIOError,
                std::format(
                    "Failed to read number of columns from metadata file: {}",
                    snp_stats.string())});
    }
    meta_stream.read(reinterpret_cast<char*>(&num_mono_snp), sizeof(int64_t));
    if (!meta_stream)
    {
        return std::unexpected(
            Error{
                ErrorCode::FileIOError,
                std::format(
                    "Failed to read number of monomorphic SNPs from metadata "
                    "file: {}",
                    snp_stats.string())});
    }

    if (rows <= 0 || cols <= 0)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidData,
                std::format(
                    "Matrix dimensions in metadata are invalid: rows={}, "
                    "cols={}",
                    rows,
                    cols)});
    }
    std::vector<int64_t> mono_indices;
    mono_indices.resize(static_cast<size_t>(num_mono_snp));

    meta_stream.read(
        reinterpret_cast<char*>(mono_indices.data()),
        sizeof(int64_t) * static_cast<size_t>(num_mono_snp));

    std::unordered_set<int64_t> mono_set(
        mono_indices.begin(), mono_indices.end());

    Eigen::VectorXd mean = Eigen::VectorXd::Zero(cols);
    Eigen::VectorXd variance = Eigen::VectorXd::Ones(cols);
    meta_stream.read(
        reinterpret_cast<char*>(mean.data()),
        sizeof(double) * static_cast<size_t>(cols));
    meta_stream.read(
        reinterpret_cast<char*>(variance.data()),
        sizeof(double) * static_cast<size_t>(cols));

    // Create memory mapping
    mio::mmap_source mmap;
    try
    {
        mmap = mio::mmap_source(bin_file.string());
    }
    catch (const std::system_error& e)
    {
        return std::unexpected(
            Error{
                ErrorCode::FileIOError,
                std::format(
                    "Failed to memory map file: {} - {}",
                    bin_file.string(),
                    e.what())});
    }

    // Verify binary file size matches expected dimensions
    const size_t expected_size = static_cast<size_t>(rows)
                                 * static_cast<size_t>(cols) * sizeof(double);
    if (mmap.size() != expected_size)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidData,
                std::format(
                    "Binary file size does not match dimensions in metadata. "
                    "Expected: {} bytes, Actual: {} bytes",
                    expected_size,
                    mmap.size())});
    }

    // Create Eigen map
    const auto* data_ptr = reinterpret_cast<const double*>(mmap.data());

#ifdef USE_AVX512
    Eigen::Map<
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::Aligned64>
        mat(data_ptr, rows, cols);
#else
    Eigen::Map<
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::Aligned32>
        mat(data_ptr, rows, cols);
#endif

    return GenotypeMap(
        std::move(mmap),
        std::move(mat),
        std::move(mono_set),
        std::move(mean),
        std::move(variance),
        rows,
        cols);
}

}  // namespace gelex
