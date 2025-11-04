#pragma once

#include <optional>
#include <vector>

#include <Eigen/Core>

#include "gelex/types/mcmc_samples.h"

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
    PosteriorSummary() = default;

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
        : coeffs(sample.coeffs[0].rows()), variance(1)
    {
    }

    PosteriorSummary coeffs;
    PosteriorSummary variance;
};

struct BaseMarkerSummary
{
    explicit BaseMarkerSummary(const BaseMarkerSamples& samples)
        : coeffs(samples.coeffs[0].rows()),
          variance(1),
          pve(samples.coeffs[0].rows())
    {
        if (!samples.tracker.empty())  // mixture model
        {
            pip = Eigen::VectorXd::Zero(samples.tracker[0].rows());
            comp_probs = Eigen::MatrixXd::Zero(
                samples.tracker[0].rows(), samples.n_proportions);
        }

        if (!samples.mixture_proportion.empty())
        {
            mixture_proportion
                = PosteriorSummary(samples.mixture_proportion[0].rows());
        }
    }

    PosteriorSummary coeffs;
    PosteriorSummary variance;
    PosteriorSummary pve;

    PosteriorSummary mixture_proportion;
    Eigen::VectorXd pip;         // Posterior inclusion probability
    Eigen::MatrixXd comp_probs;  // Per-component posterior probabilities
};

struct AdditiveSummary : BaseMarkerSummary
{
    explicit AdditiveSummary(const AdditiveSamples& samples)
        : BaseMarkerSummary(samples)
    {
    }
};

struct DominantSummary : BaseMarkerSummary
{
    explicit DominantSummary(const DominantSamples& samples)
        : BaseMarkerSummary(samples)
    {
    }
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

    const FixedSummary* fixed() const
    {
        return fixed_ ? &fixed_.value() : nullptr;
    }
    const std::vector<RandomSummary>& random() const { return random_; }
    const AdditiveSummary* additive() const
    {
        return additive_ ? &additive_.value() : nullptr;
    }
    const DominantSummary* dominant() const
    {
        return dominant_ ? &dominant_.value() : nullptr;
    }
    const PosteriorSummary& residual() const { return residual_; }

   private:
    friend class SnpEffectsWriter;

    MCMCSamples samples_;

    std::optional<FixedSummary> fixed_;
    std::vector<RandomSummary> random_;
    std::optional<AdditiveSummary> additive_;
    std::optional<DominantSummary> dominant_;
    PosteriorSummary residual_;

    double prob_;
    double phenotype_var_;

    Eigen::VectorXd additive_variances_;
    Eigen::VectorXd additive_means_;
    Eigen::VectorXd dominant_variances_;
    Eigen::VectorXd dominant_means_;
};

}  // namespace gelex
