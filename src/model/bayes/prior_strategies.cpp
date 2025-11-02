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

auto BayesADPrior::set_dominant_effect_prior(
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

auto BayesBDPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    const size_t var_size = bayes::get_cols(effect.design_matrix);
    return detail::set_pi_mixture_prior_strategy(
        effect, prior.dominant, prior, var_size);
}

// --- BayesC / BayesCD ---
auto BayesCPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    return detail::set_pi_mixture_prior_strategy(
        effect, prior.additive, prior, 1);
}

auto BayesCDPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    return detail::set_pi_mixture_prior_strategy(
        effect, prior.dominant, prior, 1);
}

auto BayesRRPrior::set_additive_effect_prior(
    bayes::AdditiveEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    return detail::set_non_mixture_prior_strategy(
        effect, prior.additive, prior);
}

auto BayesRRDPrior::set_dominant_effect_prior(
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
    return detail::set_scale_mixture_prior_strategy(
        effect, prior.additive, prior);
}

auto BayesRDPrior::set_dominant_effect_prior(
    bayes::DominantEffect& effect,
    const PriorConfig& prior) -> std::expected<void, Error>
{
    return detail::set_scale_mixture_prior_strategy(
        effect, prior.dominant, prior);
}
}  // namespace gelex
