#pragma once

#include <expected>
#include <filesystem>

#include <gelex/mio.h>
#include <Eigen/Core>
#include <unordered_set>

#include "gelex/error.h"

namespace gelex
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

    bool is_monomorphic(Eigen::Index snp_index) const
    {
        return mono_set_.contains(snp_index);
    }
    const Eigen::VectorXd& mean() const { return mean_; }
    const Eigen::VectorXd& stddev() const { return stddev_; }
    int64_t num_mono() const { return static_cast<int64_t>(mono_set_.size()); }
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
        std::unordered_set<int64_t>&& mono_set,
        Eigen::VectorXd&& mean,
        Eigen::VectorXd&& stddev,
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
    std::unordered_set<int64_t> mono_set_;
    Eigen::VectorXd mean_;
    Eigen::VectorXd stddev_;
    int64_t rows_;
    int64_t cols_;
};

}  // namespace gelex
