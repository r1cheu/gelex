#include "gelex/data/genotype_mmap.h"

#include <algorithm>
#include <format>

#include "gelex/exception.h"

namespace gelex
{

GenotypeMap::GenotypeMap(const std::filesystem::path& bin_file)
    : mat_(nullptr, 0, 0)
{
    auto snp_stats = bin_file;
    snp_stats.replace_extension(".snpstats");

    load_metadata(snp_stats);

    std::error_code ec;
    mmap_.map(bin_file.string(), ec);
    if (ec)
    {
        throw FileIOException(
            std::format(
                "Failed to mmap file {}: {}", bin_file.string(), ec.message()));
    }

    const size_t expected_size = static_cast<size_t>(rows_)
                                 * static_cast<size_t>(cols_) * sizeof(double);
    if (mmap_.size() != expected_size)
    {
        throw InvalidDataException(
            std::format(
                "Binary file size mismatch. Expected {} bytes, got {} bytes.",
                expected_size,
                mmap_.size()));
    }

    const auto* data_ptr = reinterpret_cast<const double*>(mmap_.data());

    validate_alignment(data_ptr);

    new (&mat_) MapType(data_ptr, rows_, cols_);
}

void GenotypeMap::load_metadata(const std::filesystem::path& meta_path)
{
    std::ifstream meta_stream(meta_path, std::ios::in | std::ios::binary);
    if (!meta_stream)
    {
        throw FileNotFoundException(meta_path);
    }

    rows_ = detail::read_scalar<int64_t>(meta_stream, "rows");
    cols_ = detail::read_scalar<int64_t>(meta_stream, "cols");
    auto num_mono_snp
        = detail::read_scalar<int64_t>(meta_stream, "num_mono_snp");

    if (rows_ <= 0 || cols_ <= 0)
    {
        throw InvalidDataException(
            std::format("Invalid dimensions: {}x{}", rows_, cols_));
    }

    if (num_mono_snp > 0)
    {
        mono_indices_.resize(static_cast<size_t>(num_mono_snp));
        detail::read_binary(
            meta_stream,
            mono_indices_.data(),
            mono_indices_.size(),
            "mono indices");
        std::ranges::sort(mono_indices_);  // 确保二分查找可用
    }

    mean_.resize(cols_);
    stddev_.resize(cols_);

    detail::read_binary(
        meta_stream, mean_.data(), static_cast<size_t>(cols_), "mean values");
    detail::read_binary(
        meta_stream,
        stddev_.data(),
        static_cast<size_t>(cols_),
        "stddev values");
}

bool GenotypeMap::is_monomorphic(Eigen::Index snp_index) const noexcept
{
    return std::ranges::binary_search(mono_indices_, snp_index);
}

void GenotypeMap::validate_alignment(const void* ptr)
{
    auto addr = reinterpret_cast<std::uintptr_t>(ptr);
    if (addr % ALIGNMENT_BYTES != 0)
    {
        throw InvalidDataException(
            std::format(
                "Memory alignment error. Address {} is not aligned to {} "
                "bytes. "
                "Please check mmap offset or file generation.",
                addr,
                ALIGNMENT_BYTES));
    }
}

}  // namespace gelex
