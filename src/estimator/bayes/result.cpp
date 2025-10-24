#include "gelex/estimator/bayes/result.h"

#include <optional>
#include <ranges>

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
      residual_(1),
      prob_(prob),
      phenotype_var_(model.phenotype_var())
{
    if (const auto* effect = model.additive(); effect)
    {
        additive_means_ = bayes::get_means(effect->design_matrix);
        additive_variances_ = bayes::get_variances(effect->design_matrix);
    }

    if (const auto* effect = model.dominant(); effect)
    {
        dominant_means_ = bayes::get_means(effect->design_matrix);
        dominant_variances_ = bayes::get_variances(effect->design_matrix);
    }

    if (const auto* sample = samples_.fixed(); sample)
    {
        fixed_.emplace(*sample);
    }
    for (const auto& sample : samples_.random())
    {
        random_.emplace_back(sample);
    }

    if (const auto* sample = samples_.additive(); sample)
    {
        additive_.emplace(*sample);
    }

    if (const auto* sample = samples_.dominant(); sample)
    {
        dominant_.emplace(*sample);
    }

    if (const auto* sample = samples_.pi(); sample)
    {
        pi_.emplace(*sample);
    }

    if (const auto* sample = samples_.snp_tracker(); sample)
    {
        snp_tracker_.emplace(*sample);
    }
}

void MCMCResult::compute(std::optional<double> prob)
{
    if (prob)
    {
        prob_ = prob.value();
    }

    // Use PosteriorCalculator for all computations
    if (const auto* sample = samples_.fixed(); fixed_ && sample != nullptr)
    {
        fixed_->coeffs = detail::PosteriorCalculator::compute_param_summary(
            samples_.fixed()->coeffs, prob_);
    }

    for (auto&& [result, sample] : std::views::zip(random_, samples_.random()))
    {
        result.coeffs = detail::PosteriorCalculator::compute_param_summary(
            sample.coeffs, prob_);
        result.variance = detail::PosteriorCalculator::compute_param_summary(
            sample.variance, prob_);
    }

    if (const auto* sample = samples_.additive();
        additive_ && sample != nullptr)
    {
        additive_->coeffs = detail::PosteriorCalculator::compute_param_summary(
            sample->coeffs, prob_);
        additive_->variance
            = detail::PosteriorCalculator::compute_param_summary(
                sample->variance, prob_);

        detail::PosteriorCalculator::compute_pve(
            additive_->pve,
            sample->coeffs,
            additive_variances_,
            phenotype_var_);
    }

    if (const auto* sample = samples_.dominant();
        dominant_ && sample != nullptr)
    {
        dominant_->coeffs = detail::PosteriorCalculator::compute_param_summary(
            sample->coeffs, prob_);
        dominant_->variance
            = detail::PosteriorCalculator::compute_param_summary(
                sample->variance, prob_);
        dominant_->ratios = detail::PosteriorCalculator::compute_param_summary(
            sample->ratios, prob_);

        detail::PosteriorCalculator::compute_pve(
            dominant_->pve,
            sample->coeffs,
            dominant_variances_,
            phenotype_var_);
    }

    if (const auto* sample = samples_.pi(); pi_ && sample != nullptr)
    {
        pi_->prop = detail::PosteriorCalculator::compute_param_summary(
            sample->prop, prob_);
    }

    if (const auto* sample = samples_.snp_tracker();
        snp_tracker_ && sample != nullptr)
    {
        snp_tracker_->pip
            = detail::PosteriorCalculator::compute_pip(sample->tracker);
    }

    residual_ = detail::PosteriorCalculator::compute_param_summary(
        samples_.residual().variance, prob_);
}
}  // namespace gelex
