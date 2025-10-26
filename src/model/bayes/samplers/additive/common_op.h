#pragma once

#include <Eigen/Core>
#include <cmath>

namespace gelex::detail::AdditiveSampler
{

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

inline auto compute_log_likelihood(
    double rhs,
    double marker_variance,
    double col_norm,
    double residual_variance,
    double logpi) -> double
{
    const double residual_over_marker_variance
        = residual_variance / marker_variance;
    const double percision_kernel
        = 1 / (col_norm + residual_over_marker_variance);

    const double logdetV = log((col_norm / residual_over_marker_variance) + 1);

    const double log_likelihood
        = (-0.5 * (logdetV - rhs * rhs * percision_kernel / residual_variance))
          + logpi;
    return log_likelihood;
}
}  // namespace gelex::detail::AdditiveSampler
