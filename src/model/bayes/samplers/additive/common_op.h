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
}  // namespace gelex::detail::AdditiveSampler
