#ifndef GELEX_MODEL_BAYES_SAMPLERS_MH_RRD_H_
#define GELEX_MODEL_BAYES_SAMPLERS_MH_RRD_H_

#include <random>

namespace gelex::bayes
{
struct AdditiveEffect;
struct AdditiveState;
struct DominantEffect;
struct DominantState;
struct ResidualState;
}  // namespace gelex::bayes
namespace gelex::detail::MH
{

auto RRD(
    const bayes::AdditiveEffect& add_effect,
    bayes::AdditiveState& add_state,
    const bayes::DominantEffect& dom_effect,
    bayes::DominantState& dom_state,
    bayes::ResidualState& residual,
    std::mt19937_64& rng) -> void;
}  // namespace gelex::detail::MH

#endif  // GELEX_MODEL_BAYES_SAMPLERS_MH_RRD_H_
