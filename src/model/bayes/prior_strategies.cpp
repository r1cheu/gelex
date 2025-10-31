#include "gelex/model/bayes/prior_strategies.h"

#include <expected>
#include "gelex/error.h"
#include "gelex/model/bayes/prior_constants.h"
#include "types/bayes_effects.h"

namespace gelex
{
auto BayesAPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const double target_variance
        = prior.additive.heritability * prior.phenotype_variance;
    auto init_marker_variance_result = compute_init_marker_variance(
        target_variance,
        bayes::get_matrix_ref(effect.design_matrix),
        prior_constants::NON_MIXTURE_PROPORTION);
    if (!init_marker_variance_result)
    {
        return std::unexpected(init_marker_variance_result.error());
    }
    const double init_marker_variance = *init_marker_variance_result;
    effect.init_marker_variance = init_marker_variance;
    effect.marker_variance_size = bayes::get_cols(effect.design_matrix);
    effect.marker_variance_prior
        = {prior_constants::MARKER_VARIANCE_SHAPE,
           prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
               * init_marker_variance};
    return {};
}

auto BayesADPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const double target_variance
        = prior.dominant.heritability * prior.phenotype_variance;
    auto init_marker_variance_result = compute_init_marker_variance(
        target_variance,
        bayes::get_matrix_ref(effect.design_matrix),
        prior_constants::NON_MIXTURE_PROPORTION);
    if (!init_marker_variance_result)
    {
        return std::unexpected(init_marker_variance_result.error());
    }
    const double init_marker_variance = *init_marker_variance_result;
    effect.init_marker_variance = init_marker_variance;
    effect.marker_variance_size = bayes::get_cols(effect.design_matrix);
    effect.marker_variance_prior
        = {prior_constants::MARKER_VARIANCE_SHAPE,
           prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
               * init_marker_variance};
    return {};
}

auto BayesBPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const double target_variance
        = prior.additive.heritability * prior.phenotype_variance;
    auto init_marker_variance_result = compute_init_marker_variance(
        target_variance,
        bayes::get_matrix_ref(effect.design_matrix),
        1 - prior.additive.mixture_proportions[0]);
    if (!init_marker_variance_result)
    {
        return std::unexpected(init_marker_variance_result.error());
    }
    const double init_marker_variance = *init_marker_variance_result;
    effect.init_pi.emplace(prior.additive.mixture_proportions);
    effect.init_marker_variance = init_marker_variance;
    effect.marker_variance_size = bayes::get_cols(effect.design_matrix);
    effect.marker_variance_prior
        = {prior_constants::MARKER_VARIANCE_SHAPE,
           prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
               * init_marker_variance};
    return {};
}

auto BayesBDPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const double target_variance
        = prior.dominant.heritability * prior.phenotype_variance;
    auto init_marker_variance_result = compute_init_marker_variance(
        target_variance,
        bayes::get_matrix_ref(effect.design_matrix),
        1 - prior.dominant.mixture_proportions[0]);
    if (!init_marker_variance_result)
    {
        return std::unexpected(init_marker_variance_result.error());
    }
    const double init_marker_variance = *init_marker_variance_result;
    effect.init_pi.emplace(prior.dominant.mixture_proportions);
    effect.init_marker_variance = init_marker_variance;
    effect.marker_variance_size = bayes::get_cols(effect.design_matrix);
    effect.marker_variance_prior
        = {prior_constants::MARKER_VARIANCE_SHAPE,
           prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
               * init_marker_variance};
    return {};
}

auto BayesCPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const double target_variance
        = prior.additive.heritability * prior.phenotype_variance;
    auto init_marker_variance_result = compute_init_marker_variance(
        target_variance,
        bayes::get_matrix_ref(effect.design_matrix),
        1 - prior.additive.mixture_proportions[0]);
    if (!init_marker_variance_result)
    {
        return std::unexpected(init_marker_variance_result.error());
    }
    const double init_marker_variance = *init_marker_variance_result;
    effect.init_pi.emplace(prior.additive.mixture_proportions);
    effect.init_marker_variance = init_marker_variance;
    effect.marker_variance_size = 1;
    effect.marker_variance_prior
        = {prior_constants::MARKER_VARIANCE_SHAPE,
           prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
               * init_marker_variance};
    return {};
}

auto BayesCDPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const double target_variance
        = prior.dominant.heritability * prior.phenotype_variance;
    auto init_marker_variance_result = compute_init_marker_variance(
        target_variance,
        bayes::get_matrix_ref(effect.design_matrix),
        1 - prior.dominant.mixture_proportions[0]);
    if (!init_marker_variance_result)
    {
        return std::unexpected(init_marker_variance_result.error());
    }
    const double init_marker_variance = *init_marker_variance_result;
    effect.init_pi.emplace(prior.dominant.mixture_proportions);
    effect.init_marker_variance = init_marker_variance;
    effect.marker_variance_size = 1;
    effect.marker_variance_prior
        = {prior_constants::MARKER_VARIANCE_SHAPE,
           prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
               * init_marker_variance};
    return {};
}

auto BayesRRPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const double target_variance
        = prior.additive.heritability * prior.phenotype_variance;
    auto init_marker_variance_result = compute_init_marker_variance(
        target_variance,
        bayes::get_matrix_ref(effect.design_matrix),
        prior_constants::NON_MIXTURE_PROPORTION);
    if (!init_marker_variance_result)
    {
        return std::unexpected(init_marker_variance_result.error());
    }
    const double init_marker_variance = *init_marker_variance_result;
    effect.init_marker_variance = init_marker_variance;
    effect.marker_variance_size = bayes::get_cols(effect.design_matrix);
    effect.marker_variance_prior
        = {prior_constants::MARKER_VARIANCE_SHAPE,
           prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
               * init_marker_variance};
    return {};
}

auto BayesRRDPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const double target_variance
        = prior.dominant.heritability * prior.phenotype_variance;
    auto init_marker_variance_result = compute_init_marker_variance(
        target_variance,
        bayes::get_matrix_ref(effect.design_matrix),
        prior_constants::NON_MIXTURE_PROPORTION);
    if (!init_marker_variance_result)
    {
        return std::unexpected(init_marker_variance_result.error());
    }
    const double init_marker_variance = *init_marker_variance_result;
    effect.init_marker_variance = init_marker_variance;
    effect.marker_variance_size = bayes::get_cols(effect.design_matrix);
    effect.marker_variance_prior
        = {prior_constants::MARKER_VARIANCE_SHAPE,
           prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
               * init_marker_variance};
    return {};
}

auto BayesRPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const double target_variance
        = prior.additive.heritability * prior.phenotype_variance;
    auto init_marker_variance_result = compute_init_marker_variance(
        target_variance,
        bayes::get_matrix_ref(effect.design_matrix),
        1 - prior.additive.mixture_proportions[0]);
    if (!init_marker_variance_result)
    {
        return std::unexpected(init_marker_variance_result.error());
    }
    const double init_marker_variance = *init_marker_variance_result;
    effect.init_marker_variance = init_marker_variance;
    effect.scale.emplace(prior.additive.mixture_scales);
    effect.marker_variance_size = 1;
    effect.marker_variance_prior
        = {prior_constants::MARKER_VARIANCE_SHAPE,
           prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
               * init_marker_variance};
    return {};
}

auto BayesRDPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const double target_variance
        = prior.dominant.heritability * prior.phenotype_variance;
    auto init_marker_variance_result = compute_init_marker_variance(
        target_variance,
        bayes::get_matrix_ref(effect.design_matrix),
        1 - prior.dominant.mixture_proportions[0]);
    if (!init_marker_variance_result)
    {
        return std::unexpected(init_marker_variance_result.error());
    }
    const double init_marker_variance = *init_marker_variance_result;
    effect.init_marker_variance = init_marker_variance;
    effect.scale.emplace(prior.dominant.mixture_scales);
    effect.marker_variance_size = 1;
    effect.marker_variance_prior
        = {prior_constants::MARKER_VARIANCE_SHAPE,
           prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
               * init_marker_variance};
    return {};
}
}  // namespace gelex
