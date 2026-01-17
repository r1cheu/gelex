#ifndef GELEX_MODEL_BAYES_SAMPLERS_DOMINANT_H_
#define GELEX_MODEL_BAYES_SAMPLERS_DOMINANT_H_

#include <random>

namespace gelex
{
class BayesModel;
class BayesState;
}  // namespace gelex

namespace gelex::detail::DominantSampler
{

struct A
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

struct B
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

struct C
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

struct R
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

struct RR
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

struct Coeff
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

struct RatioMean
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

struct RatioVar
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

}  // namespace gelex::detail::DominantSampler

#endif  // GELEX_MODEL_BAYES_SAMPLERS_DOMINANT_H_
