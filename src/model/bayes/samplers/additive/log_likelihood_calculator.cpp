#include "log_likelihood_calculator.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace gelex::detail::AdditiveSampler
{

LogLikelihoodCalculator::LogLikelihoodCalculator(
    double col_norm,
    double rhs,
    double residual_variance,
    const Eigen::VectorXd& logpi)
    : col_norm_(col_norm),
      rhs_(rhs),
      residual_variance_(residual_variance),
      logpi_(logpi)
{
    const auto num_components = logpi.size();
    marker_variances_.reserve(num_components);
    precision_kernels_.reserve(num_components);
    log_likelihoods_.reserve(num_components);
}

void LogLikelihoodCalculator::add_component(double marker_variance)
{
    marker_variances_.push_back(marker_variance);
    precision_kernels_.push_back(compute_precision_kernel(marker_variance));
}

void LogLikelihoodCalculator::compute_probabilities()
{
    const int num_components = static_cast<int>(marker_variances_.size());

    // Compute log-likelihoods for all components
    log_likelihoods_.resize(num_components);
    for (int k = 0; k < num_components; ++k)
    {
        log_likelihoods_[k] = compute_log_likelihood(k);
    }

    // Convert to probabilities using log-sum-exp trick for numerical stability
    probabilities_.resize(num_components);

    // Find maximum log-likelihood for numerical stability
    double max_log_likelihood = std::ranges::max(log_likelihoods_);

    // Compute unnormalized probabilities
    double sum_exp = 0.0;
    for (int k = 0; k < num_components; ++k)
    {
        probabilities_(k) = std::exp(log_likelihoods_[k] - max_log_likelihood);
        sum_exp += probabilities_(k);
    }

    // Normalize probabilities
    probabilities_ /= sum_exp;
}

double LogLikelihoodCalculator::precision_kernel(int component_index) const
{
    return precision_kernels_[component_index];
}

double LogLikelihoodCalculator::posterior_mean(int component_index) const
{
    return rhs_ * precision_kernels_[component_index];
}

double LogLikelihoodCalculator::posterior_stddev(int component_index) const
{
    return std::sqrt(residual_variance_ * precision_kernels_[component_index]);
}

double LogLikelihoodCalculator::compute_precision_kernel(
    double marker_variance) const
{
    const double residual_over_marker_variance
        = residual_variance_ / marker_variance;
    return 1.0 / (col_norm_ + residual_over_marker_variance);
}

double LogLikelihoodCalculator::compute_logdetV(double marker_variance) const
{
    return std::log((marker_variance * col_norm_ / residual_variance_) + 1.0);
}

double LogLikelihoodCalculator::compute_log_likelihood(
    int component_index) const
{
    if (component_index == 0)
    {
        // Component 0 corresponds to exclusion (coefficient = 0)
        return logpi_(0);
    }
    // For non-zero components, compute full log-likelihood
    const double precision_kernel = precision_kernels_[component_index];
    const double logdetV = compute_logdetV(marker_variances_[component_index]);

    return (-0.5
            * (logdetV - (rhs_ * rhs_ * precision_kernel / residual_variance_)))
           + logpi_(component_index);
}

}  // namespace gelex::detail::AdditiveSampler
