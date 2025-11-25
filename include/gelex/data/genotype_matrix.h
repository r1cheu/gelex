#pragma once

#include <algorithm>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

class GenotypeMatrix
{
   public:
    GenotypeMatrix(
        Eigen::MatrixXd&& data,
        std::vector<int64_t>&& mono_indices,
        Eigen::VectorXd&& mean,
        Eigen::VectorXd&& stddev);

    GenotypeMatrix(const GenotypeMatrix&) = delete;
    GenotypeMatrix(GenotypeMatrix&&) noexcept = default;
    GenotypeMatrix& operator=(const GenotypeMatrix&) = delete;
    GenotypeMatrix& operator=(GenotypeMatrix&&) noexcept = default;
    ~GenotypeMatrix() = default;

    [[nodiscard]] const Eigen::MatrixXd& matrix() const noexcept
    {
        return data_;
    }

    [[nodiscard]] bool is_monomorphic(Eigen::Index marker_idx) const noexcept
    {
        return std::ranges::binary_search(mono_indices_, marker_idx);
    }

    [[nodiscard]] const Eigen::VectorXd& mean() const noexcept { return mean_; }
    [[nodiscard]] const Eigen::VectorXd& stddev() const noexcept
    {
        return stddev_;
    }

    [[nodiscard]] int64_t num_mono() const noexcept
    {
        return static_cast<int64_t>(mono_indices_.size());
    }
    [[nodiscard]] int64_t rows() const noexcept { return data_.rows(); }
    [[nodiscard]] int64_t cols() const noexcept { return data_.cols(); }

   private:
    Eigen::MatrixXd data_;
    std::vector<int64_t> mono_indices_;
    Eigen::VectorXd mean_;
    Eigen::VectorXd stddev_;

    void validate_dimensions() const;
};

}  // namespace gelex
