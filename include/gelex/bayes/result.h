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
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/detail/result_factory.h"
#include "gelex/bayes/draws.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
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

template <GeneticModeSet Modes>
class VarianceSummaryResult
{
   public:
    VarianceSummaryResult(
        HomogeneousModeValues<Modes, ScalarPosteriorResult> explained_variance,
        HomogeneousModeValues<Modes, ScalarPosteriorResult> heritability,
        ScalarPosteriorResult total_explained_variance,
        ScalarPosteriorResult total_heritability)
        : explained_variance_{std::move(explained_variance)},
          heritability_{std::move(heritability)},
          total_explained_variance_{std::move(total_explained_variance)},
          total_heritability_{std::move(total_heritability)}
    {
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto explained_variance() const noexcept
        -> const ScalarPosteriorResult&
    {
        return explained_variance_.template get<Mode>();
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto heritability() const noexcept
        -> const ScalarPosteriorResult&
    {
        return heritability_.template get<Mode>();
    }

    [[nodiscard]] auto explained_variances() const noexcept
        -> const HomogeneousModeValues<Modes, ScalarPosteriorResult>&
    {
        return explained_variance_;
    }

    [[nodiscard]] auto heritabilities() const noexcept
        -> const HomogeneousModeValues<Modes, ScalarPosteriorResult>&
    {
        return heritability_;
    }

    [[nodiscard]] auto total_explained_variance() const noexcept
        -> const ScalarPosteriorResult&
    {
        return total_explained_variance_;
    }

    [[nodiscard]] auto total_heritability() const noexcept
        -> const ScalarPosteriorResult&
    {
        return total_heritability_;
    }

   private:
    HomogeneousModeValues<Modes, ScalarPosteriorResult> explained_variance_;
    HomogeneousModeValues<Modes, ScalarPosteriorResult> heritability_;
    ScalarPosteriorResult total_explained_variance_;
    ScalarPosteriorResult total_heritability_;
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
    using genetic_result_type = detail::genetic_result_t<GeneticPrior>;

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

    [[nodiscard]] auto genetic() const noexcept -> const genetic_result_type&
    {
        return genetic_;
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
        genetic_result_type genetic,
        ScalarPosteriorResult residual,
        VarianceSummaryResult<modes> variance_summary)
        : fixed_{std::move(fixed)},
          random_{std::move(random)},
          genetic_{std::move(genetic)},
          residual_{std::move(residual)},
          variance_summary_{std::move(variance_summary)}
    {
    }

    template <typename T>
    friend auto make_result(const BayesModel& model, const BayesDraws<T>& draws)
        -> BayesResult<T>;

    CoefficientPosteriorResult fixed_;
    std::vector<RandomEffectResult> random_;
    genetic_result_type genetic_;
    ScalarPosteriorResult residual_;
    VarianceSummaryResult<modes> variance_summary_;
};

namespace detail
{

template <GeneticModeSet Modes>
auto make_variance_summary_result(const VarianceSummaryDraws<Modes>& draws)
    -> VarianceSummaryResult<Modes>
{
    auto explained_variance = generate_mode_values<Modes>(
        [&]<GeneticMode Mode>()
        {
            return make_scalar_result(
                draws.template explained_variance<Mode>());
        });
    auto heritability = generate_mode_values<Modes>(
        [&]<GeneticMode Mode>()
        { return make_scalar_result(draws.template heritability<Mode>()); });
    auto total_explained_variance
        = make_scalar_result(draws.total_explained_variance());
    auto total_heritability = make_scalar_result(draws.total_heritability());
    return VarianceSummaryResult<Modes>{
        std::move(explained_variance),
        std::move(heritability),
        std::move(total_explained_variance),
        std::move(total_heritability)};
}

}  // namespace detail

template <typename GeneticPrior>
auto make_result(const BayesModel& model, const BayesDraws<GeneticPrior>& draws)
    -> BayesResult<GeneticPrior>
{
    auto fixed = detail::make_coefficient_result(
        draws.fixed(), model.fixed().column_names(), model.fixed().X().cols());

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

    std::vector<RandomEffectResult> random;
    random.reserve(random_designs.size());
    for (const auto [index, design] : random_designs | std::views::enumerate)
    {
        const auto& block = random_draws[static_cast<std::size_t>(index)];
        random.emplace_back(
            detail::make_coefficient_result(
                block.coefficients(), design.column_names(), design.X().cols()),
            detail::make_scalar_result(block.variance()),
            detail::make_scalar_result(block.explained_variance()));
    }

    auto genetic = detail::make_genetic_result(
        std::type_identity<GeneticPrior>{},
        draws.genetic(),
        model.genetic().cols());
    auto residual = detail::make_scalar_result(draws.residual());
    auto variance_summary
        = detail::make_variance_summary_result(draws.variance_summary());

    return BayesResult<GeneticPrior>{
        std::move(fixed),
        std::move(random),
        std::move(genetic),
        std::move(residual),
        std::move(variance_summary)};
}

}  // namespace gelex

#endif  // GELEX_BAYES_RESULT_H_
