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

#ifndef GELEX_BAYES_STATE_H_
#define GELEX_BAYES_STATE_H_

#include <Eigen/Core>
#include <fmt/format.h>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include "gelex/bayes/detail/state_factory.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/exception.h"

namespace gelex
{

struct FixedEffectState
{
    Eigen::VectorXd coefficients;
};

struct RandomEffectState
{
    Eigen::VectorXd coefficients;
    Eigen::VectorXd fitted_values;
    double variance{};
};

struct ResidualState
{
    Eigen::VectorXd adjusted_response;
    double variance{};
};

template <typename GeneticPrior>
class BayesState;

template <typename GeneticPrior>
[[nodiscard]] auto make_state(
    const BayesPrior<GeneticPrior>& prior,
    const BayesModel& model) -> BayesState<GeneticPrior>;

template <typename GeneticPrior>
class BayesState
{
   public:
    using genetic_prior_type = GeneticPrior;
    using genetic_state_type = detail::genetic_state_t<GeneticPrior>;

    [[nodiscard]] auto fixed() noexcept -> FixedEffectState& { return fixed_; }
    [[nodiscard]] auto fixed() const noexcept -> const FixedEffectState&
    {
        return fixed_;
    }

    [[nodiscard]] auto random() noexcept -> std::span<RandomEffectState>
    {
        return random_;
    }
    [[nodiscard]] auto random() const noexcept
        -> std::span<const RandomEffectState>
    {
        return random_;
    }

    [[nodiscard]] auto genetic() noexcept -> genetic_state_type&
    {
        return genetic_;
    }
    [[nodiscard]] auto genetic() const noexcept -> const genetic_state_type&
    {
        return genetic_;
    }

    [[nodiscard]] auto residual() noexcept -> ResidualState&
    {
        return residual_;
    }
    [[nodiscard]] auto residual() const noexcept -> const ResidualState&
    {
        return residual_;
    }

   private:
    BayesState(
        FixedEffectState fixed,
        std::vector<RandomEffectState> random,
        genetic_state_type genetic,
        ResidualState residual)
        : fixed_{std::move(fixed)},
          random_{std::move(random)},
          genetic_{std::move(genetic)},
          residual_{std::move(residual)}
    {
    }

    template <typename T>
    friend auto make_state(const BayesPrior<T>& prior, const BayesModel& model)
        -> BayesState<T>;

    FixedEffectState fixed_;
    std::vector<RandomEffectState> random_;
    genetic_state_type genetic_;
    ResidualState residual_;
};

template <typename GeneticPrior>
[[nodiscard]] auto make_state(
    const BayesPrior<GeneticPrior>& prior,
    const BayesModel& model) -> BayesState<GeneticPrior>
{
    const auto random_prior = prior.random();
    const auto random_design = model.random();
    if (random_prior.size() != random_design.size())
    {
        throw GelexException(
            fmt::format(
                "random prior/design count mismatch: {} != {}",
                random_prior.size(),
                random_design.size()));
    }

    std::vector<RandomEffectState> random;
    random.reserve(random_design.size());
    for (const auto& [parameter, design] :
         std::views::zip(random_prior, random_design))
    {
        random.push_back(
            RandomEffectState{
                .coefficients = Eigen::VectorXd::Zero(design.X().cols()),
                .fitted_values = Eigen::VectorXd::Zero(design.X().rows()),
                .variance = parameter.initial});
    }

    return BayesState<GeneticPrior>{
        FixedEffectState{
            .coefficients = Eigen::VectorXd::Zero(model.fixed().X().cols())},
        std::move(random),
        detail::make_state(prior.genetic(), model.genetic()),
        ResidualState{
            .adjusted_response = model.phenotype(),
            .variance = prior.residual().initial}};
}

}  // namespace gelex

#endif  // GELEX_BAYES_STATE_H_
