#pragma once

#include <expected>
#include <unordered_set>

#include <Eigen/Core>

#include "gelex/error.h"

namespace gelex
{

/**
 * @class GenotypeMatrix
 * @brief In-memory genotype matrix storage (alternative to memory-mapped
 * GenotypeMap)
 *
 * Stores genotype data directly in memory using Eigen::MatrixXd.
 * Provides the same interface as GenotypeMap for drop-in replacement.
 * Suitable for smaller datasets that fit in RAM or when file I/O should be
 * minimized.
 */
class GenotypeMatrix
{
   public:
    /**
     * @brief Create from existing matrix and metadata
     *
     * @param data Genotype matrix (rows=samples, cols=markers)
     * @param mono_set Set of monomorphic marker indices
     * @param mean Mean values for each marker
     * @param variance Variance values for each marker
     * @return GenotypeMatrix instance
     */
    static auto create(
        Eigen::MatrixXd&& data,
        std::unordered_set<int64_t>&& mono_set,
        Eigen::VectorXd&& mean,
        Eigen::VectorXd&& stddev) -> GenotypeMatrix;

    /**
     * @brief Get reference to the genotype matrix
     */
    const Eigen::MatrixXd& matrix() const { return data_; }

    /**
     * @brief Check if a marker is monomorphic
     */
    bool is_monomorphic(Eigen::Index marker_idx) const
    {
        return mono_set_.contains(marker_idx);
    }

    /**
     * @brief Get mean values for all markers
     */
    const Eigen::VectorXd& mean() const { return mean_; }

    /**
     * @brief Get variance values for all markers
     */
    const Eigen::VectorXd& stddev() const { return stddev_; }

    /**
     * @brief Get number of monomorphic markers
     */
    int64_t num_mono() const { return static_cast<int64_t>(mono_set_.size()); }

    /**
     * @brief Get number of samples (rows)
     */
    int64_t rows() const { return data_.rows(); }

    /**
     * @brief Get number of markers (columns)
     */
    int64_t cols() const { return data_.cols(); }

   private:
    explicit GenotypeMatrix(
        Eigen::MatrixXd&& data,
        std::unordered_set<int64_t>&& mono_set,
        Eigen::VectorXd&& mean,
        Eigen::VectorXd&& stddev)
        : data_(std::move(data)),
          mono_set_(std::move(mono_set)),
          mean_(std::move(mean)),
          stddev_(std::move(stddev))
    {
    }

    Eigen::MatrixXd data_;
    std::unordered_set<int64_t> mono_set_;
    Eigen::VectorXd mean_;
    Eigen::VectorXd stddev_;
};

}  // namespace gelex
