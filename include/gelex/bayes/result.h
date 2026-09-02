/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GELEX_BAYES_RESULT_H_
#define GELEX_BAYES_RESULT_H_

#include <cstddef>
#include <fmt/format.h>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/detail/result_factory.h"
#include "gelex/bayes/draws.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/variance/result.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

class RandomEffectResult
{
   public:
    RandomEffectResult(
        CoefficientPosteriorResult coefficients,
        ScalarPosteriorResult variance,
        ScalarPosteriorResult explained_variance)
        : coefficients_{std::move(coefficients)},
          variance_{std::move(variance)},
          explained_variance_{std::move(explained_variance)}
    {
    }

    [[nodiscard]] auto coefficients() const noexcept
        -> const CoefficientPosteriorResult&
    {
        return coefficients_;
    }

    [[nodiscard]] auto variance() const noexcept -> const ScalarPosteriorResult&
    {
        return variance_;
    }

    [[nodiscard]] auto explained_variance() const noexcept
        -> const ScalarPosteriorResult&
    {
        return explained_variance_;
    }

   private:
    CoefficientPosteriorResult coefficients_;
    ScalarPosteriorResult variance_;
    ScalarPosteriorResult explained_variance_;
};

template <typename GeneticPrior>
class BayesResult;

template <typename GeneticPrior>
[[nodiscard]] auto make_result(
    const BayesModel& model,
    const BayesDraws<GeneticPrior>& draws) -> BayesResult<GeneticPrior>;

template <typename GeneticPrior>
class BayesResult
{
   public:
    using genetic_prior_type = GeneticPrior;
    using genetic_parameter_result_type
        = decltype(detail::make_genetic_parameters(
            std::declval<const typename BayesDraws<
                GeneticPrior>::genetic_draws_type&>()));
    using marker_effect_result_type = decltype(detail::make_marker_effects(
        std::declval<const BayesModel&>(),
        std::declval<
            const typename BayesDraws<GeneticPrior>::genetic_draws_type&>()));

    static constexpr GeneticModeSet modes = GeneticPrior::modes;

    [[nodiscard]] auto fixed() const noexcept
        -> const CoefficientPosteriorResult&
    {
        return fixed_;
    }

    [[nodiscard]] auto random() const noexcept
        -> std::span<const RandomEffectResult>
    {
        return random_;
    }

    [[nodiscard]] auto genetic_parameters() const noexcept
        -> const genetic_parameter_result_type&
    {
        return genetic_parameters_;
    }

    [[nodiscard]] auto marker_effects() const noexcept
        -> const marker_effect_result_type&
    {
        return marker_effects_;
    }

    [[nodiscard]] auto residual() const noexcept -> const ScalarPosteriorResult&
    {
        return residual_;
    }

    [[nodiscard]] auto variance_summary() const noexcept
        -> const VarianceSummaryResult<modes>&
    {
        return variance_summary_;
    }

   private:
    BayesResult(
        CoefficientPosteriorResult fixed,
        std::vector<RandomEffectResult> random,
        genetic_parameter_result_type genetic_parameters,
        marker_effect_result_type marker_effects,
        ScalarPosteriorResult residual,
        VarianceSummaryResult<modes> variance_summary)
        : fixed_{std::move(fixed)},
          random_{std::move(random)},
          genetic_parameters_{std::move(genetic_parameters)},
          marker_effects_{std::move(marker_effects)},
          residual_{std::move(residual)},
          variance_summary_{std::move(variance_summary)}
    {
    }

    template <typename T>
    friend auto make_result(const BayesModel& model, const BayesDraws<T>& draws)
        -> BayesResult<T>;

    CoefficientPosteriorResult fixed_;
    std::vector<RandomEffectResult> random_;
    genetic_parameter_result_type genetic_parameters_;
    marker_effect_result_type marker_effects_;
    ScalarPosteriorResult residual_;
    VarianceSummaryResult<modes> variance_summary_;
};

template <typename GeneticPrior>
auto make_result(const BayesModel& model, const BayesDraws<GeneticPrior>& draws)
    -> BayesResult<GeneticPrior>
{
    auto fixed
        = detail::make_result(draws.fixed(), model.fixed().column_names());

    const auto random_designs = model.random();
    const auto random_draws = draws.random();
    if (random_designs.size() != random_draws.size())
    {
        throw GelexException(
            fmt::format(
                "model has {} random blocks but draws have {}",
                random_designs.size(),
                random_draws.size()));
    }
    const auto random_explained_variance = draws.variance_summary().random();

    std::vector<RandomEffectResult> random;
    random.reserve(random_designs.size());
    for (const auto [index, design] : random_designs | std::views::enumerate)
    {
        const auto& block = random_draws[static_cast<std::size_t>(index)];
        const auto& explained_variance
            = random_explained_variance[static_cast<std::size_t>(index)];
        random.emplace_back(
            detail::make_result(block.coefficients(), design.column_names()),
            detail::make_result(block.variance()),
            detail::make_result(explained_variance));
    }

    auto genetic_parameters = detail::make_genetic_parameters(draws.genetic());
    auto marker_effects = detail::make_marker_effects(model, draws.genetic());
    auto residual = detail::make_result(draws.residual());
    auto variance_summary = detail::make_result(draws.variance_summary());

    return BayesResult<GeneticPrior>{
        std::move(fixed),
        std::move(random),
        std::move(genetic_parameters),
        std::move(marker_effects),
        std::move(residual),
        std::move(variance_summary)};
}

}  // namespace gelex

#endif  // GELEX_BAYES_RESULT_H_
