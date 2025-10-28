#pragma once

#include <cassert>

#include <mkl.h>

#include <Eigen/Core>
#include <cmath>

namespace gelex::detail::AdditiveSampler
{

struct LikelihoodParams
{
    double log_likelihood{0.0};
    double precision_kernel{0.0};
    double residual_over_marker_variance{0.0};
};

struct PosteriorParams
{
    double mean{0.0};
    double stddev{0.0};
    double log_likelihood_kernel{0.0};
};

template <typename DerivedX, typename DerivedY>
inline double mkl_ddot(
    const Eigen::DenseBase<DerivedX>& x,
    const Eigen::DenseBase<DerivedY>& y)
{
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedX);
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedY);
    assert(x.size() == y.size() && "mkl_ddot: vector sizes do not match.");

    const MKL_INT n = static_cast<MKL_INT>(x.size());
    const MKL_INT incx = 1;
    const MKL_INT incy = 1;

    return ddot(&n, x.derived().data(), &incx, y.derived().data(), &incy);
}

template <typename DerivedX, typename DerivedY>
inline void mkl_daxpy(
    double alpha,
    const Eigen::DenseBase<DerivedX>& x,
    Eigen::DenseBase<DerivedY>& y)
{
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedX);
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedY);
    assert(x.size() == y.size() && "mkl_daxpy: vector sizes do not match.");

    const MKL_INT n = static_cast<MKL_INT>(x.size());
    const double a = alpha;
    const MKL_INT incx = 1;
    const MKL_INT incy = 1;

    daxpy(&n, &a, x.derived().data(), &incx, y.derived().data(), &incy);
}

inline auto update_residual_and_gebv(
    Eigen::Ref<Eigen::VectorXd> y_adj,
    Eigen::Ref<Eigen::VectorXd> gebv,
    const Eigen::Ref<const Eigen::VectorXd>& col,
    double old_value,
    double new_value = 0.0) -> void
{
    const double diff = old_value - new_value;
    if (fabs(diff) > std::numeric_limits<double>::epsilon())
    {
        mkl_daxpy(diff, col, y_adj);
        mkl_daxpy(-diff, col, gebv);
    }
}

inline auto compute_likelihood_params(
    double rhs,
    double marker_variance,
    double col_norm,
    double residual_variance,
    double logpi) -> LikelihoodParams
{
    const double res_over_marker_var = residual_variance / marker_variance;
    const double precision_k = 1.0 / (col_norm + res_over_marker_var);

    const double logdetV = std::log((col_norm / res_over_marker_var) + 1.0);

    const double log_like
        = (-0.5 * (logdetV - rhs * rhs * precision_k / residual_variance))
          + logpi;

    return {log_like, precision_k, res_over_marker_var};
}

inline auto compute_posterior_params_core(
    double rhs,
    double col_norm,
    double residual_variance,
    double res_over_marker_var) -> PosteriorParams
{
    const double precision_kernel = 1.0 / (col_norm + res_over_marker_var);

    const double post_mean = rhs * precision_kernel;
    const double post_stddev = std::sqrt(residual_variance * precision_kernel);

    const double logdetV = std::log(col_norm / res_over_marker_var + 1.0);

    const double log_like_kernel
        = -0.5 * (logdetV - post_mean * rhs / residual_variance);

    return {post_mean, post_stddev, log_like_kernel};
}

inline auto compute_posterior_params(
    double rhs,
    double marker_variance_i,
    double col_norm,
    double residual_variance) -> PosteriorParams
{
    const double res_over_marker_var = residual_variance / marker_variance_i;
    return compute_posterior_params_core(
        rhs, col_norm, residual_variance, res_over_marker_var);
}
}  // namespace gelex::detail::AdditiveSampler
