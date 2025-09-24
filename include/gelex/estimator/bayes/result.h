#pragma once
#include <filesystem>
#include <memory>
#include <optional>
#include <vector>

#include <Eigen/Core>

#include "gelex/estimator/bayes/samples.h"

namespace gelex
{
struct PosteriorSummary
{
    explicit PosteriorSummary(Eigen::Index n_params)
        : mean(Eigen::VectorXd::Zero(n_params)),
          stddev(Eigen::VectorXd::Zero(n_params)),
          hpdi_low(Eigen::VectorXd::Zero(n_params)),
          hpdi_high(Eigen::VectorXd::Zero(n_params)),
          ess(Eigen::VectorXd::Zero(n_params)),
          rhat(Eigen::VectorXd::Zero(n_params))
    {
    }

    Eigen::Index size() const { return mean.size(); }

    Eigen::VectorXd mean;
    Eigen::VectorXd stddev;
    Eigen::VectorXd hpdi_low;
    Eigen::VectorXd hpdi_high;
    Eigen::VectorXd ess;   // Effective sample size
    Eigen::VectorXd rhat;  // Gelman-Rubin diagnostic

    std::array<double, 6> operator[](Eigen::Index idx) const
    {
        return std::array<double, 6>{
            mean(idx),
            stddev(idx),
            hpdi_low(idx),
            hpdi_high(idx),
            ess(idx),
            rhat(idx)};
    }
};

struct RandomSummary
{
    RandomSummary(Eigen::Index n_coeff, Eigen::Index n_sigma)
        : coeff(n_coeff), sigma(n_sigma)
    {
    }
    explicit RandomSummary(const RandomSamples& sample)
        : RandomSummary(sample.coeffs[0].rows(), sample.sigmas[0].rows())
    {
    }

    PosteriorSummary coeff;
    PosteriorSummary sigma;
};

struct AdditiveSummary : RandomSummary
{
    AdditiveSummary(
        Eigen::Index n_coeff,
        Eigen::Index n_sigma,
        Eigen::Index n_variance)
        : RandomSummary(n_coeff, n_sigma), variance(n_variance)
    {
    }

    explicit AdditiveSummary(const AdditiveSamples& samples)
        : AdditiveSummary(
              samples.coeffs[0].rows(),
              samples.sigmas[0].rows(),
              samples.variance[0].rows())
    {
    }

    PosteriorSummary variance;
};

struct DominantSummary
{
    explicit DominantSummary(Eigen::Index n_coeff, Eigen::Index n_variance)
        : coeff(n_coeff), variance(n_variance)
    {
    }

    explicit DominantSummary(const DominantSamples& samples)
        : DominantSummary(samples.coeffs[0].rows(), samples.variance[0].rows())
    {
    }

    PosteriorSummary coeff;
    PosteriorSummary variance;
};

class MCMCResult
{
   public:
    explicit MCMCResult(MCMCSamples&& samples, double prob = 0.9);

    /**
     * @brief Compute posterior statistics.
     *
     * If prob is provided, uses it as the probability threshold for
     * computation. Otherwise, uses default prob.
     *
     * @param prob Optional probability threshold for computation.
     */
    void compute(std::optional<double> prob = std::nullopt);
    void save(const std::filesystem::path& prefix) const;

    const PosteriorSummary& fixed() const { return *fixed_; }
    const std::vector<RandomSummary>& random() const { return random_; }
    const std::unique_ptr<AdditiveSummary>& additive() const
    {
        return additive_;
    }
    const std::unique_ptr<DominantSummary>& dominant() const
    {
        return dominant_;
    }
    const PosteriorSummary& residual() const { return residual_; }

   private:
    /**
     * @brief Compute posterior summary statistics (mean, std, HPDI, ess, rhat).
     *
     * @param summary Destination PosteriorSummary to populate.
     * @param samples  MCMC samples organized by chain (rows=params,
     * cols=draws).
     * @param prob     Probability mass for the HPDI (e.g. 0.95).
     */
    static void compute_summary_statistics(
        PosteriorSummary& summary,
        const Samples& samples,
        double prob);

    /**
     * @brief Compute basic posterior summary statistics (mean, std).
     *
     * @param summary Destination PosteriorSummary to populate.
     * @param samples  MCMC samples organized by chain (rows=params,
     * cols=draws).
     */
    static void compute_summary_statistics(
        PosteriorSummary& summary,
        const Samples& samples);

    MCMCSamples samples_;

    Eigen::VectorXd add_eff_;
    Eigen::VectorXd dom_eff_;

    std::unique_ptr<PosteriorSummary> fixed_;
    std::vector<RandomSummary> random_;
    std::unique_ptr<AdditiveSummary> additive_ = nullptr;
    std::unique_ptr<DominantSummary> dominant_ = nullptr;
    PosteriorSummary residual_;

    double prob_;
};

}  // namespace gelex
