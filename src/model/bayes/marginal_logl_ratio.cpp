#include "../src/model/bayes/marginal_logl_ratio.h"

#include <cmath>

#include <Eigen/Core>

namespace gelex
{
namespace bayes
{

double MarginalLogLRatioCalculator::calculate_log_likelihood(
    double M_j_dot_y_star,
    double M_j_squared_norm,
    double alpha_j,
    double sigma_alpha_k,
    double sigma_e,
    double pi_k)
{
    // Calculate the two terms in the log-likelihood formula
    const double term1
        = std::log((M_j_squared_norm * sigma_alpha_k / sigma_e) + 1.0);

    const double numerator
        = std::pow(M_j_dot_y_star + M_j_squared_norm * alpha_j, 2);
    const double denominator
        = (M_j_squared_norm + sigma_e / sigma_alpha_k) * sigma_e;
    const double term2 = numerator / denominator;

    // Combine terms and add log(pi_k)
    const double log_likelihood = -0.5 * (term1 - term2) + std::log(pi_k);

    return log_likelihood;
}

double MarginalLogLRatioCalculator::calculate_distribution_probability(
    const Eigen::Ref<const Eigen::VectorXd>& log_likelihoods,
    int k)
{
    if (k < 0 || k >= log_likelihoods.size())
    {
        return 0.0;
    }

    double sum = 0.0;
    const double L_k = log_likelihoods(k);

    for (int i = 0; i < log_likelihoods.size(); ++i)
    {
        sum += std::exp(log_likelihoods(i) - L_k);
    }

    return 1.0 / sum;
}

Eigen::VectorXd MarginalLogLRatioCalculator::calculate_all_probabilities(
    const Eigen::Ref<const Eigen::VectorXd>& log_likelihoods)
{
    const int n_distributions = log_likelihoods.size();
    Eigen::VectorXd probabilities(n_distributions);

    for (int k = 0; k < n_distributions; ++k)
    {
        probabilities(k)
            = calculate_distribution_probability(log_likelihoods, k);
    }

    return probabilities;
}

Eigen::VectorXd MarginalLogLRatioCalculator::calculate_all_log_likelihoods(
    double M_j_dot_y_star,
    double M_j_squared_norm,
    double alpha_j,
    const Eigen::Ref<const Eigen::VectorXd>& sigma_alpha,
    double sigma_e,
    const Eigen::Ref<const Eigen::VectorXd>& pi)
{
    const int n_distributions = sigma_alpha.size();
    Eigen::VectorXd log_likelihoods(n_distributions);

    for (int k = 0; k < n_distributions; ++k)
    {
        log_likelihoods(k) = calculate_log_likelihood(
            M_j_dot_y_star,
            M_j_squared_norm,
            alpha_j,
            sigma_alpha(k),
            sigma_e,
            pi(k));
    }

    return log_likelihoods;
}

}  // namespace bayes
}  // namespace gelex
