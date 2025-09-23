#pragma once

#include <expected>
#include <filesystem>

#include <gelex/mio.h>
#include <Eigen/Core>

#include "gelex/error.h"

namespace gelex::detail
{

#ifdef USE_AVX512
static constexpr size_t ALIGNMENT = 64;
#else
static constexpr size_t ALIGNMENT = 32;
#endif

class GenotypeMap
{
   public:
    static auto create(const std::filesystem::path& bin_file)
        -> std::expected<GenotypeMap, Error>;

    const Eigen::Map<
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
#ifdef USE_AVX512
        Eigen::Aligned64
#else
        Eigen::Aligned32
#endif
        >&
    matrix() const
    {
        return mat_;
    }

    int64_t rows() const { return rows_; }
    int64_t cols() const { return cols_; }

   private:
    explicit GenotypeMap(
        mio::mmap_source&& mmap,
#ifdef USE_AVX512
        Eigen::Map<
            const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
            Eigen::Aligned64>&& mat,
#else
        Eigen::Map<
            const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
            Eigen::Aligned32>&& mat,
#endif
        int64_t rows,
        int64_t cols)
        : mmap_(std::move(mmap)), mat_(std::move(mat)), rows_(rows), cols_(cols)
    {
    }

    mio::mmap_source mmap_;
#ifdef USE_AVX512
    Eigen::Map<
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::Aligned64>
        mat_;
#else
    Eigen::Map<
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::Aligned32>
        mat_;
#endif
    int64_t rows_;
    int64_t cols_;
};

}  // namespace gelex::detail
