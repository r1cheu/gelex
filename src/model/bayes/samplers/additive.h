#ifndef GELEX_MODEL_BAYES_SAMPLERS_ADDITIVE_H_
#define GELEX_MODEL_BAYES_SAMPLERS_ADDITIVE_H_

#include <random>

namespace gelex
{
class BayesModel;
class BayesState;
}  // namespace gelex

namespace gelex::detail::AdditiveSampler
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

struct RRD
{
    auto operator()(
        const BayesModel& model,
        BayesState& states,
        std::mt19937_64& rng) const -> void;
};

}  // namespace gelex::detail::AdditiveSampler

#endif  // GELEX_MODEL_BAYES_SAMPLERS_ADDITIVE_H_
