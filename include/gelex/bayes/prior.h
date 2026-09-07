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

#include <cmath>
#include <fmt/format.h>
#include <span>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/construction.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/exception.h"
#include "gelex/infra/var.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

template <typename GeneticPrior>
class BayesPrior
{
   public:
    using genetic_prior_type = GeneticPrior;
    BayesPrior(
        std::vector<VarianceParameter> random,
        GeneticPrior genetic,
        VarianceParameter residual)
        : random_{std::move(random)},
          genetic_{std::move(genetic)},
          residual_{residual}
    {
    }

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
    std::vector<VarianceParameter> random_;
    GeneticPrior genetic_;
    VarianceParameter residual_;
};

GELEX_NAMESPACE_BEGIN(detail)
inline auto random_projection_variance(const bayes::RandomDesign& design)
    -> double
{
    const double variance = matvar(design.X(), VarNormType::Population).sum();
    if (!std::isfinite(variance) || variance <= 0.0)
    {
        throw GelexException(
            fmt::format(
                "random design '{}' projection variance must be finite and "
                "positive, got {}",
                design.name(),
                variance));
    }
    return variance;
}

inline auto make_random_variance_parameters(
    const BayesModel& model,
    const VarianceBudget& budget,
    double phenotype_variance) -> std::vector<VarianceParameter>
{
    const auto designs = model.random();
    const double share = budget.random();
    if (designs.empty())
    {
        if (share != 0.0)
        {
            throw GelexException(
                "random variance share must be zero when the model has no "
                "random designs");
        }
        return {};
    }
    if (share <= 0.0)
    {
        throw GelexException(
            "random variance share must be positive when the model has "
            "random designs");
    }

    const double block_target
        = phenotype_variance * share / static_cast<double>(designs.size());
    std::vector<VarianceParameter> parameters;
    parameters.reserve(designs.size());
    for (const auto& design : designs)
    {
        const double initial
            = block_target / random_projection_variance(design);
        parameters.push_back(make_mean_calibrated_variance_parameter(initial));
    }
    return parameters;
}
GELEX_NAMESPACE_END(detail)

template <GeneticModeSet Modes, typename GeneticSpec>
[[nodiscard]] auto make_prior(
    const BayesRecipe<Modes, GeneticSpec>& recipe,
    const BayesModel& model)
{
    const double phenotype_variance = model.phenotype_variance();
    const detail::MarkerVarianceCalibrator calibrator{model, recipe.variance()};
    auto genetic = detail::make_prior(recipe.genetic_spec(), calibrator);
    auto random = detail::make_random_variance_parameters(
        model, recipe.variance(), phenotype_variance);
    auto residual = detail::make_mean_calibrated_variance_parameter(
        phenotype_variance * recipe.variance().residual());
    return BayesPrior<decltype(genetic)>{
        std::move(random), std::move(genetic), residual};
}

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_PRIOR_H_
