#pragma once

#include <cmath>

#include <Eigen/Core>

namespace gelex
{
namespace bayes
{

class MarginalLogLRatioCalculator
{
   public:
    MarginalLogLRatioCalculator() = default;

    /**
     * @brief Calculate the log-likelihood for a genetic marker being in
     * distribution k
     *
     * @param M_j_dot_y_star Dot product M_j^T y*
     * @param M_j_squared_norm Squared norm M_j^T M_j (col_norm)
     * @param alpha_j Current effect size for marker j
     * @param sigma_alpha_k Variance parameter for distribution k
     * @param sigma_e Residual variance
     * @param pi_k Probability for distribution k
     * @return double Log-likelihood L_π_k
     */
    static double calculate_log_likelihood(
        double M_j_dot_y_star,
        double col_norm,
        double coeff,
        double sigma,
        double sigma_e,
        double pi_k);

    /**
     * @brief Calculate probability that marker j is in distribution k
     *
     * @param log_likelihoods Vector of log-likelihoods for all distributions
     * @param k Index of the distribution
     * @return double Probability P(π_k)
     */
    static double calculate_distribution_probability(
        const Eigen::Ref<const Eigen::VectorXd>& log_likelihoods,
        int k);

    /**
     * @brief Calculate probabilities for all distributions
     *
     * @param log_likelihoods Vector of log-likelihoods for all distributions
     * @return Eigen::VectorXd Vector of probabilities for each distribution
     */
    static Eigen::VectorXd calculate_all_probabilities(
        const Eigen::Ref<const Eigen::VectorXd>& log_likelihoods);

    /**
     * @brief Calculate all log-likelihoods for a genetic marker across
     * distributions
     *
     * @param M_j_dot_y_star Dot product M_j^T y*
     * @param M_j_squared_norm Squared norm M_j^T M_j (col_norm)
     * @param alpha_j Current effect size for marker j
     * @param sigma_alpha Vector of variance parameters for each distribution
     * @param sigma_e Residual variance
     * @param pi Vector of probabilities for each distribution
     * @return Eigen::VectorXd Vector of log-likelihoods for each distribution
     */
    static Eigen::VectorXd calculate_all_log_likelihoods(
        double M_j_dot_y_star,
        double M_j_squared_norm,
        double alpha_j,
        const Eigen::Ref<const Eigen::VectorXd>& sigma_alpha,
        double sigma_e,
        const Eigen::Ref<const Eigen::VectorXd>& pi);
};

}  // namespace bayes
}  // namespace gelex
