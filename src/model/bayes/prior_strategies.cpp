#include "gelex/model/bayes/prior_strategies.h"

#include "../src/types/bayes_effects.h"
#include "gelex/exception.h"
#include "gelex/model/bayes/prior_strategy.h"

namespace gelex
{
auto BayesAPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> void
{
    detail::set_non_mixture_prior_strategy(effect, prior.additive, prior);
}

auto BayesAdPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> void
{
    detail::set_non_mixture_prior_strategy(effect, prior.dominant, prior);
}

auto BayesBPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> void
{
    const size_t var_size = bayes::get_cols(effect.design_matrix);
    detail::set_pi_mixture_prior_strategy(
        effect, prior.additive, prior, var_size);
}

auto BayesBdPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> void
{
    const size_t var_size = bayes::get_cols(effect.design_matrix);
    detail::set_pi_mixture_prior_strategy(
        effect, prior.dominant, prior, var_size);
}

auto BayesBpiPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> void
{
    const size_t var_size = bayes::get_cols(effect.design_matrix);
    detail::set_pi_mixture_prior_strategy(
        effect, prior.additive, prior, var_size);
    effect.estimate_pi = true;
}

auto BayesBdpiPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> void
{
    const size_t var_size = bayes::get_cols(effect.design_matrix);
    detail::set_pi_mixture_prior_strategy(
        effect, prior.dominant, prior, var_size);
    effect.estimate_pi = true;
}

auto BayesCPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> void
{
    detail::set_pi_mixture_prior_strategy(effect, prior.additive, prior, 1);
}

auto BayesCdPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> void
{
    detail::set_pi_mixture_prior_strategy(effect, prior.dominant, prior, 1);
}

auto BayesCpiPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> void
{
    detail::set_pi_mixture_prior_strategy(effect, prior.additive, prior, 1);
    effect.estimate_pi = true;
}

auto BayesCdpiPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> void
{
    detail::set_pi_mixture_prior_strategy(effect, prior.dominant, prior, 1);
    effect.estimate_pi = true;
}

auto BayesRRPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> void
{
    detail::set_non_mixture_prior_strategy(effect, prior.additive, prior);
}

auto BayesRRdPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> void
{
    detail::set_non_mixture_prior_strategy(effect, prior.dominant, prior);
}

// --- BayesR / BayesRD ---
auto BayesRPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> void
{
    detail::set_scale_mixture_prior_strategy(effect, prior.additive, prior);
    effect.estimate_pi = true;
}

auto BayesRdPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> void
{
    detail::set_scale_mixture_prior_strategy(effect, prior.dominant, prior);
    effect.estimate_pi = true;
}
}  // namespace gelex
