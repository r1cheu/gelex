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
        : FixedSummary(sample.coeffs[0].rows())
    {
    }
    PosteriorSummary coeffs;

   protected:
    explicit FixedSummary(Eigen::Index n_coeff) : coeffs(n_coeff) {}
};

struct RandomSummary : FixedSummary
{
    explicit RandomSummary(const RandomSamples& sample)
        : RandomSummary(sample.coeffs[0].rows(), 1)
    {
    }

    PosteriorSummary variance;

   protected:
    RandomSummary(Eigen::Index n_coeff, Eigen::Index n_variance)
        : FixedSummary(n_coeff), variance(n_variance)
    {
    }
};

struct AdditiveSummary : RandomSummary
{
    explicit AdditiveSummary(const AdditiveSamples& samples)
        : AdditiveSummary(samples.coeffs[0].rows())
    {
    }
    PosteriorSummary pve;

   protected:
    explicit AdditiveSummary(Eigen::Index n_coeff)
        : RandomSummary(n_coeff, 1), pve(n_coeff)
    {
    }
};

struct DominantSummary : AdditiveSummary
{
    explicit DominantSummary(const DominantSamples& samples)
        : DominantSummary(samples.coeffs[0].rows())
    {
    }
    PosteriorSummary ratios;

   protected:
    explicit DominantSummary(Eigen::Index n_coeff)
        : AdditiveSummary(n_coeff), ratios(n_coeff)
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
    friend class MCMCResultWriter;

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
