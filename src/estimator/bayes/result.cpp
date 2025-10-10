#include "gelex/estimator/bayes/result.h"

#include <cstddef>
#include <memory>
#include <optional>

#include <Eigen/Core>

#include "../src/estimator/bayes/posterior_calculator.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{

using Eigen::Index;
using Eigen::VectorXd;

MCMCResult::MCMCResult(
    MCMCSamples&& samples,
    const BayesModel& model,
    double prob)
    : samples_(std::move(samples)),
      prob_(prob),
      phenotype_var_(model.phenotype_var()),
      residual_(1)
{
    // Extract phenotype variance from model

    // Extract additive variances and means if additive effect exists
    if (model.additive())
    {
        const auto& additive_matrix = model.additive()->design_matrix;
        additive_variances_ = additive_matrix.variance();
        additive_means_ = additive_matrix.mean();
    }

    // Extract dominant variances and means if dominant effect exists
    if (model.dominant())
    {
        const auto& dominant_matrix = model.dominant()->design_matrix;
        dominant_variances_ = dominant_matrix.variance();
        dominant_means_ = dominant_matrix.mean();
    }

    fixed_ = samples_.fixed() ? std::make_unique<FixedSummary>(samples_.fixed())
                              : nullptr;
    additive_ = samples_.additive()
                    ? std::make_unique<AdditiveSummary>(samples_.additive())
                    : nullptr;
    dominant_ = samples_.dominant()
                    ? std::make_unique<DominantSummary>(samples_.dominant())
                    : nullptr;

    for (const auto& sample : samples_.random().chain_samples)
    {
        random_.emplace_back(sample);
    }
}

void MCMCResult::compute(std::optional<double> prob)
{
    if (prob)
    {
        prob_ = prob.value();
    }

    // Use PosteriorCalculator for all computations
    if (samples_.fixed())
    {
        fixed_->coeffs = detail::PosteriorCalculator::compute_param_summary(
            samples_.fixed().coeffs, prob_);
    }

    for (size_t i = 0; i < samples_.random().size(); ++i)
    {
        random_[i].coeffs = detail::PosteriorCalculator::compute_param_summary(
            samples_.random().chain_samples[i].coeffs, prob_);
        random_[i].effect_variance
            = detail::PosteriorCalculator::compute_param_summary(
                samples_.random().chain_samples[i].effect_variance, prob_);
    }

    if (samples_.additive())
    {
        additive_->coeffs = detail::PosteriorCalculator::compute_snp_summary(
            samples_.additive().coeffs);
        additive_->effect_variance
            = detail::PosteriorCalculator::compute_param_summary(
                samples_.additive().effect_variance, prob_);

        detail::PosteriorCalculator::compute_pve(
            additive_->pve,
            samples_.additive().coeffs,
            additive_variances_,
            phenotype_var_);
    }

    if (samples_.dominant())
    {
        dominant_->coeffs = detail::PosteriorCalculator::compute_snp_summary(
            samples_.dominant().coeffs);
        dominant_->effect_variance
            = detail::PosteriorCalculator::compute_param_summary(
                samples_.dominant().effect_variance, prob_);
        detail::PosteriorCalculator::compute_pve(
            dominant_->pve,
            samples_.dominant().coeffs,
            dominant_variances_,
            phenotype_var_);
    }

    residual_ = detail::PosteriorCalculator::compute_param_summary(
        samples_.residual().variance, prob_);
}
}  // namespace gelex
