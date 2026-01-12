#include "gelex/optim/convergence_checker.h"

#include <cmath>
#include <limits>

namespace gelex
{
bool ConvergenceChecker::is_converged(
    const Eigen::Ref<const Eigen::VectorXd>& new_sigma,
    double new_loglike)
{
    double sigma_diff = compute_sigma_diff(new_sigma);
    double loglike_diff = compute_loglike_diff(new_loglike);

    bool negative = loglike_diff < 0.0;
    loglike_diff = std::abs(loglike_diff);

    if (sigma_diff < tol_
        && (loglike_diff < 1e-4 || (negative && loglike_diff < 1e-2)))
    {
        converged_ = true;
    }
    update(new_sigma, new_loglike);

    return converged_;
}

double ConvergenceChecker::compute_sigma_diff(
    const Eigen::Ref<const Eigen::VectorXd>& new_sigma) const
{
    if (old_sigma_.size() == 0 || old_sigma_.size() != new_sigma.size())
    {
        return std::numeric_limits<double>::infinity();
    }
    return (new_sigma - old_sigma_).norm() / new_sigma.norm();
}

double ConvergenceChecker::compute_loglike_diff(double new_loglike) const
{
    return new_loglike - old_loglike_;
}
}  // namespace gelex
