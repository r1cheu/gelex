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

#include "gelex/simulate/genetic_value_scaler.h"

#include <cmath>

#include <Eigen/Core>

#include "gelex/infra/stats/descriptive.h"
#include "gelex/types/sim_types.h"

namespace gelex
{

auto GeneticValueScaler::scale(
    const GeneticValues& additive_values,
    GeneticValues* dominance_values,
    Eigen::Ref<Eigen::VectorXd> residual) const -> void
{
    const double additive_variance = detail::var(additive_values.gebv)(0);
    if (additive_variance <= 0.0 || h2_ <= 0.0)
    {
        return;
    }

    const double d2 = d2_.value_or(0.0);

    if (dominance_values != nullptr && d2 > 0.0)
    {
        const double dominance_raw_var = detail::var(dominance_values->gebv)(0);
        if (dominance_raw_var > 0.0)
        {
            const double target_dom_var = additive_variance * d2 / h2_;
            const double dom_scale
                = std::sqrt(target_dom_var / dominance_raw_var);
            dominance_values->gebv *= dom_scale;
            dominance_values->coeff *= dom_scale;
        }
    }

    const double residual_variance = additive_variance * (1.0 - h2_ - d2) / h2_;
    if (residual_variance > 0.0)
    {
        residual *= std::sqrt(residual_variance);
    }
    else
    {
        residual.setZero();
    }
}

}  // namespace gelex
