#include "gelex/model/bayes/prior_strategies.h"

#include <expected>

#include "gelex/error.h"
#include "gelex/model/bayes/prior_strategy.h"
#include "types/bayes_effects.h"

namespace gelex
{
auto BayesAPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    return detail::set_non_mixture_prior_strategy(
        effect, prior.additive, prior);
}

auto BayesAdPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    return detail::set_non_mixture_prior_strategy(
        effect, prior.dominant, prior);
}

auto BayesBPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const size_t var_size = bayes::get_cols(effect.design_matrix);
    return detail::set_pi_mixture_prior_strategy(
        effect, prior.additive, prior, var_size);
}

auto BayesBdPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const size_t var_size = bayes::get_cols(effect.design_matrix);
    return detail::set_pi_mixture_prior_strategy(
        effect, prior.dominant, prior, var_size);
}

auto BayesBpiPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const size_t var_size = bayes::get_cols(effect.design_matrix);
    auto result = detail::set_pi_mixture_prior_strategy(
        effect, prior.additive, prior, var_size);
    if (result)
    {
        effect.estimate_pi = true;
    }
    return result;
}

auto BayesBdpiPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const size_t var_size = bayes::get_cols(effect.design_matrix);
    auto result = detail::set_pi_mixture_prior_strategy(
        effect, prior.dominant, prior, var_size);
    if (result)
    {
        effect.estimate_pi = true;
    }
    return result;
}

// --- BayesC / BayesCD ---
auto BayesCPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    return detail::set_pi_mixture_prior_strategy(
        effect, prior.additive, prior, 1);
}

auto BayesCdPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    return detail::set_pi_mixture_prior_strategy(
        effect, prior.dominant, prior, 1);
}

auto BayesCpiPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    auto result = detail::set_pi_mixture_prior_strategy(
        effect, prior.additive, prior, 1);
    if (result)
    {
        effect.estimate_pi = true;
    }
    return result;
}

auto BayesCdpiPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    auto result = detail::set_pi_mixture_prior_strategy(
        effect, prior.dominant, prior, 1);
    if (result)
    {
        effect.estimate_pi = true;
    }
    return result;
}

auto BayesRRPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    return detail::set_non_mixture_prior_strategy(
        effect, prior.additive, prior);
}

auto BayesRRdPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    return detail::set_non_mixture_prior_strategy(
        effect, prior.dominant, prior);
}

// --- BayesR / BayesRD ---
auto BayesRPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    auto result = detail::set_scale_mixture_prior_strategy(
        effect, prior.additive, prior);
    if (result)
    {
        effect.estimate_pi = true;
    }
    return result;
}

auto BayesRdPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    auto result = detail::set_scale_mixture_prior_strategy(
        effect, prior.dominant, prior);
    if (result)
    {
        effect.estimate_pi = true;
    }
    return result;
}
}  // namespace gelex
