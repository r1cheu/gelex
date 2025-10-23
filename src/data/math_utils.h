#pragma once

#include <cmath>    // for erf, sqrt
#include <numbers>  // for std::numbers::sqrt2

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace gelex
{
namespace detail
{
Eigen::RowVectorXd centralize(Eigen::Ref<Eigen::MatrixXd> x);
std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd> standardize(
    Eigen::Ref<Eigen::MatrixXd> x);
Eigen::VectorXd sum_square(const Eigen::Ref<const Eigen::MatrixXd>& mat);
Eigen::VectorXd sum_square(const Eigen::Ref<Eigen::SparseMatrix<double>>& mat);
Eigen::VectorXd cols_var(const Eigen::Ref<const Eigen::MatrixXd>& mat);

// Normal CDF for any mean (mu) and stddev (sigma)
inline double normal_cdf(double x, double mu = 0.0, double sigma = 1.0)
{
    return 0.5 * (1.0 + std::erf((x - mu) / (sigma * std::numbers::sqrt2)));
}

template <typename Derived>
Eigen::VectorXd var(
    const Eigen::DenseBase<Derived>& a,
    Eigen::Index norm_type = 1,
    Eigen::Index axis = 0)
{
    const Eigen::Index n = (axis == 0) ? a.cols() : a.rows();
    const Eigen::Index ddof = (norm_type == 0) ? 0 : 1;

    Eigen::VectorXd result(n);

    if (axis == 0)
    {
#pragma omp parallel for default(none) shared(n, a, result, ddof)
        for (Eigen::Index i = 0; i < n; ++i)
        {
            auto col = a.col(i);
            double mean_val = col.mean();
            double sum_sq = (col.array() - mean_val).square().sum();
            result(i) = sum_sq / static_cast<double>(col.size() - ddof);
        }
    }
    else
    {
        // Compute variance along rows (axis=1)
#pragma omp parallel for default(none) shared(n, a, result, ddof)
        for (Eigen::Index i = 0; i < n; ++i)
        {
            auto row = a.row(i);
            double mean_val = row.mean();
            double sum_sq = (row.array() - mean_val).square().sum();
            result(i) = sum_sq / static_cast<double>(row.size() - ddof);
        }
    }

    return result;
}

}  // namespace detail
}  // namespace gelex
