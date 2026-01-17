#ifndef GELEX_MODEL_BAYES_SAMPLERS_GIBBS_GIBBS_CONCEPT_H_
#define GELEX_MODEL_BAYES_SAMPLERS_GIBBS_GIBBS_CONCEPT_H_

#include <concepts>

namespace gelex::bayes
{
struct AdditiveEffect;
struct AdditiveState;
struct DominantEffect;
struct DominantState;
}  // namespace gelex::bayes

namespace gelex::detail::Gibbs
{

template <typename E, typename S>
concept IsValidEffectStatePair
    = (std::same_as<std::remove_cvref_t<E>, bayes::AdditiveEffect>
       && std::same_as<std::remove_cvref_t<S>, bayes::AdditiveState>)
      || (std::same_as<std::remove_cvref_t<E>, bayes::DominantEffect>
          && std::same_as<std::remove_cvref_t<S>, bayes::DominantState>);

}  // namespace gelex::detail::Gibbs

#endif  // GELEX_MODEL_BAYES_SAMPLERS_GIBBS_GIBBS_CONCEPT_H_
