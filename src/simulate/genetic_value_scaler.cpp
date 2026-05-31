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

#include "gelex/infra/stats/detail/var.h"
#include "gelex/simulate/sim_types.h"

namespace gelex
{

auto GeneticValueScaler::scale(
    GeneticValues* additive_values,
    GeneticValues* dominance_values,
    Eigen::Ref<Eigen::VectorXd> residual) const -> void
{
    const bool has_add = additive_values != nullptr && h2_ && *h2_ > 0.0;
    const bool has_dom = dominance_values != nullptr && d2_ && *d2_ > 0.0;

    double ref_var = 0.0;
    double trait_frac = 0.0;
    if (has_add)
    {
        ref_var = stats::detail::vecvar(additive_values->gebv);
        trait_frac = *h2_;
    }
    else if (has_dom)
    {
        ref_var = stats::detail::vecvar(dominance_values->gebv);
        trait_frac = *d2_;
    }
    else
    {
        return;
    }

    if (ref_var <= 0.0)
    {
        return;
    }

    const double total_var = ref_var / trait_frac;

    if (has_add && has_dom)
    {
        const double dominance_raw_var
            = stats::detail::vecvar(dominance_values->gebv);
        if (dominance_raw_var > 0.0)
        {
            const double target_dom_var = total_var * *d2_;
            const double dom_scale
                = std::sqrt(target_dom_var / dominance_raw_var);
            dominance_values->gebv *= dom_scale;
            dominance_values->coeff *= dom_scale;
        }
    }

    const double h2 = has_add ? *h2_ : 0.0;
    const double d2 = has_dom ? *d2_ : 0.0;
    const double residual_variance = total_var * (1.0 - h2 - d2);
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
