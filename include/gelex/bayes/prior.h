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

#ifndef GELEX_BAYES_PRIOR_H_
#define GELEX_BAYES_PRIOR_H_

#include <span>
#include <utility>
#include <vector>

#include "gelex/bayes/detail/prior_factory.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/variance/parameter.h"

namespace gelex
{

template <typename GeneticPrior>
class BayesPrior;

template <GeneticModeSet Modes, typename GeneticFamily>
[[nodiscard]] auto make_prior(
    const BayesRecipe<Modes, GeneticFamily>& recipe,
    const BayesModel& model);

template <typename GeneticPrior>
class BayesPrior
{
   public:
    using genetic_prior_type = GeneticPrior;

    [[nodiscard]] auto random() const noexcept
        -> std::span<const VarianceParameter>
    {
        return random_;
    }

    [[nodiscard]] auto genetic() const noexcept -> const GeneticPrior&
    {
        return genetic_;
    }

    [[nodiscard]] auto residual() const noexcept -> const VarianceParameter&
    {
        return residual_;
    }

   private:
    BayesPrior(
        std::vector<VarianceParameter> random,
        GeneticPrior genetic,
        VarianceParameter residual)
        : random_{std::move(random)},
          genetic_{std::move(genetic)},
          residual_{residual}
    {
    }

    template <GeneticModeSet Modes, typename GeneticFamily>
    friend auto make_prior(
        const BayesRecipe<Modes, GeneticFamily>& recipe,
        const BayesModel& model);

    std::vector<VarianceParameter> random_;
    GeneticPrior genetic_;
    VarianceParameter residual_;
};

template <GeneticModeSet Modes, typename GeneticFamily>
[[nodiscard]] auto make_prior(
    const BayesRecipe<Modes, GeneticFamily>& recipe,
    const BayesModel& model)
{
    const double phenotype_variance = model.phenotype_variance();
    const detail::MarkerVarianceCalibrator calibrator{model, recipe.variance()};
    auto genetic = detail::make_prior<Modes>(
        GeneticFamily{}, recipe.genetic_spec(), calibrator);
    auto random = detail::make_random_variance_parameters(
        model, recipe.variance(), phenotype_variance);
    auto residual = detail::make_mean_calibrated_variance_parameter(
        phenotype_variance * recipe.variance().residual());
    return BayesPrior<decltype(genetic)>{
        std::move(random), std::move(genetic), residual};
}

}  // namespace gelex

#endif  // GELEX_BAYES_PRIOR_H_
