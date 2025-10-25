#pragma once

#include <Eigen/Core>

namespace gelex::detail::AdditiveSampler
{

/**
 * @class LogLikelihoodCalculator
 * @brief Calculates log-likelihoods and component probabilities for Bayesian
 * mixture models
 *
 * This class abstracts the complex log-likelihood calculations used in Bayesian
 * mixture models (like BayesC and BayesR) to avoid code duplication and provide
 * efficient computation of component probabilities and posterior parameters.
 */
class LogLikelihoodCalculator
{
   public:
    /**
     * @brief Construct a new Log Likelihood Calculator object
     *
     * @param col_norm The squared norm of the design matrix column
     * @param rhs The right-hand side value (col.dot(y_adj) + col_norm *
     * old_value)
     * @param residual_variance The residual variance
     * @param logpi Log of mixture proportions for each component
     */
    LogLikelihoodCalculator(
        double col_norm,
        double rhs,
        double residual_variance,
        const Eigen::VectorXd& logpi);

    /**
     * @brief Add a component with specific marker variance
     *
     * @param marker_variance The marker variance for this component
     */
    void add_component(double marker_variance);

    /**
     * @brief get the posterior parameters based accept_probility
     *
     * @param marker_variance The marker variance for this component
     */
    std::pair<double, double> get_posterior_params(double accept_prob) const;
    /**
     * @brief Compute all component probabilities
     *
     * This method calculates the log-likelihoods for all components and
     * converts them to normalized probabilities suitable for sampling.
     */
    void compute_probabilities();

    /**
     * @brief Get the component probabilities
     *
     * @return const Eigen::VectorXd& Vector of normalized probabilities for
     * each component
     */
    const Eigen::VectorXd& probabilities() const { return probabilities_; }

    /**
     * @brief Get the precision kernel for a specific component
     *
     * @param component_index The index of the component
     * @return double The precision kernel value
     */
    double precision_kernel(int component_index) const;

    /**
     * @brief Get the posterior mean for a specific component
     *
     * @param component_index The index of the component
     * @return double The posterior mean value
     */
    double posterior_mean(int component_index) const;

    /**
     * @brief Get the posterior standard deviation for a specific component
     *
     * @param component_index The index of the component
     * @return double The posterior standard deviation
     */
    double posterior_stddev(int component_index) const;

    /**
     * @brief Get the number of components
     *
     * @return int Number of components added
     */
    int num_components() const
    {
        return static_cast<int>(marker_variances_.size());
    }

   private:
    double col_norm_;
    double rhs_;
    double residual_variance_;
    Eigen::VectorXd logpi_;

    std::vector<double> marker_variances_;
    std::vector<double> precision_kernels_;
    std::vector<double> log_likelihoods_;
    Eigen::VectorXd probabilities_;

    /**
     * @brief Compute precision kernel for a given marker variance
     *
     * @param marker_variance The marker variance
     * @return double The precision kernel value
     */
    double compute_precision_kernel(double marker_variance) const;

    /**
     * @brief Compute log determinant for a given marker variance
     *
     * @param marker_variance The marker variance
     * @return double The log determinant value
     */
    double compute_logdetV(double marker_variance) const;

    /**
     * @brief Compute log-likelihood for a component
     *
     * @param component_index The index of the component
     * @return double The log-likelihood value
     */
    double compute_log_likelihood(int component_index) const;
};

}  // namespace gelex::detail::AdditiveSampler
