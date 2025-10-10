#pragma once

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
};

struct FixedSummary
{
    explicit FixedSummary(const FixedSamples& sample)
        : coeffs(sample.coeffs[0].rows())
    {
    }
    PosteriorSummary coeffs;
};

struct RandomSummary
{
    explicit RandomSummary(const RandomSamples& sample)
        : RandomSummary(
              sample.coeffs[0].rows(),
              sample.effect_variance[0].rows())
    {
    }

    PosteriorSummary coeffs;
    PosteriorSummary effect_variance;

   protected:
    RandomSummary(Eigen::Index n_coeff, Eigen::Index n_variance)
        : coeffs(n_coeff), effect_variance(n_variance)
    {
    }
};

struct AdditiveSummary : RandomSummary
{
    explicit AdditiveSummary(const AdditiveSamples& samples)
        : AdditiveSummary(
              samples.coeffs[0].rows(),
              samples.effect_variance[0].rows(),
              samples.marker_variance[0].rows())
    {
    }
    PosteriorSummary marker_variance;
    PosteriorSummary pve;

   protected:
    AdditiveSummary(
        Eigen::Index n_coeff,
        Eigen::Index n_effect_variance,
        Eigen::Index n_marker_variance)
        : RandomSummary(n_coeff, n_effect_variance),
          marker_variance(n_marker_variance),
          pve(n_coeff)
    {
    }
};

struct DominantSummary : RandomSummary
{
    explicit DominantSummary(const DominantSamples& samples)
        : RandomSummary(
              samples.coeffs[0].rows(),
              samples.effect_variance[0].rows()),
          pve(samples.coeffs[0].rows())
    {
    }

    PosteriorSummary pve;
};

class MCMCResult
{
   public:
    explicit MCMCResult(
        MCMCSamples&& samples,
        const BayesModel& model,
        double prob = 0.9);

    /**
     * @brief Compute posterior statistics.
     *
     * If prob is provided, uses it as the probability threshold for
     * computation. Otherwise, uses default prob.
     *
     * @param prob Optional probability threshold for computation.
     */
    void compute(std::optional<double> prob = std::nullopt);

    const std::unique_ptr<FixedSummary>& fixed() const { return fixed_; }
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
    friend class MCMCResultWriter;

    MCMCSamples samples_;

    std::unique_ptr<FixedSummary> fixed_ = nullptr;
    std::vector<RandomSummary> random_;
    std::unique_ptr<AdditiveSummary> additive_ = nullptr;
    std::unique_ptr<DominantSummary> dominant_ = nullptr;
    PosteriorSummary residual_;

    double prob_;
    double phenotype_var_;
    Eigen::VectorXd additive_variances_;
    Eigen::VectorXd additive_means_;
    Eigen::VectorXd dominant_variances_;
    Eigen::VectorXd dominant_means_;
};

}  // namespace gelex
