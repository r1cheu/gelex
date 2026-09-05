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

#ifndef GELEX_BAYES_DETAIL_PRIOR_FACTORY_H_
#define GELEX_BAYES_DETAIL_PRIOR_FACTORY_H_

#include <cmath>
#include <fmt/format.h>
#include <vector>

#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/bayes/variance/detail/calibration.h"
#include "gelex/exception.h"
#include "gelex/infra/var.h"

namespace gelex
{

namespace detail
{

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

}  // namespace detail

}  // namespace gelex

#endif  // GELEX_BAYES_DETAIL_PRIOR_FACTORY_H_
