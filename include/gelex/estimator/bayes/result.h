#pragma once

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

    explicit FixedSummary(Eigen::Index n_coeff) : coeffs(n_coeff) {}

    PosteriorSummary coeffs;
};

struct RandomSummary
{
    explicit RandomSummary(const RandomSamples& sample)
        : coeffs(sample.coeffs[0].rows()), variance(1)
    {
    }

    RandomSummary(Eigen::Index n_coeff, Eigen::Index n_variance)
        : coeffs(n_coeff), variance(n_variance)
    {
    }

    PosteriorSummary coeffs;
    PosteriorSummary variance;
};

struct AdditiveSummary
{
    explicit AdditiveSummary(const AdditiveSamples& samples)
        : coeffs(samples.coeffs[0].rows()),
          variance(1),
          pve(samples.coeffs[0].rows())
    {
    }

    explicit AdditiveSummary(Eigen::Index n_coeff)
        : coeffs(n_coeff), variance(1), pve(n_coeff)
    {
    }

    PosteriorSummary coeffs;
    PosteriorSummary variance;
    PosteriorSummary pve;
};

struct PiSummary
{
    explicit PiSummary(const PiSamples& samples) : prop(samples.prop[0].rows())
    {
    }

    explicit PiSummary(Eigen::Index n_props) : prop(n_props) {}

    PosteriorSummary prop;
};

struct SnpTrackerSummary
{
    explicit SnpTrackerSummary(const SnpTrackerSamples& samples)
        : pip(Eigen::VectorXd::Zero(samples.tracker[0].rows()))
    {
    }

    explicit SnpTrackerSummary(Eigen::Index n_snps)
        : pip(Eigen::VectorXd::Zero(n_snps))
    {
    }

    Eigen::VectorXd pip;         // Posterior inclusion probability
    Eigen::MatrixXd comp_probs;  // Per-component posterior probabilities
                                 // (n_snps x n_components)
};

struct DominantSummary
{
    explicit DominantSummary(const DominantSamples& samples)
        : coeffs(samples.coeffs[0].rows()),
          variance(1),
          pve(samples.coeffs[0].rows()),
          ratios(samples.coeffs[0].rows())
    {
    }

    explicit DominantSummary(Eigen::Index n_coeff)
        : coeffs(n_coeff), variance(1), pve(n_coeff), ratios(n_coeff)
    {
    }

    PosteriorSummary coeffs;
    PosteriorSummary variance;
    PosteriorSummary pve;
    PosteriorSummary ratios;
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
    const PiSummary* pi() const { return pi_ ? &pi_.value() : nullptr; }
    const SnpTrackerSummary* snp_tracker() const
    {
        return snp_tracker_ ? &snp_tracker_.value() : nullptr;
    }

   private:
    friend class MCMCResultWriter;

    MCMCSamples samples_;

    std::optional<FixedSummary> fixed_;
    std::vector<RandomSummary> random_;
    std::optional<AdditiveSummary> additive_;
    std::optional<DominantSummary> dominant_;
    PosteriorSummary residual_;
    std::optional<PiSummary> pi_;
    std::optional<SnpTrackerSummary> snp_tracker_;

    double prob_;
    double phenotype_var_;

    Eigen::VectorXd additive_variances_;
    Eigen::VectorXd additive_means_;
    Eigen::VectorXd dominant_variances_;
    Eigen::VectorXd dominant_means_;
};

}  // namespace gelex
