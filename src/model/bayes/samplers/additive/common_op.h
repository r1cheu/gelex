#pragma once

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
        y_adj.noalias() += col * diff;
        gebv.noalias() -= col * diff;
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
