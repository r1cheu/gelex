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

#include "gelex/optim/constrain.h"

#include <Eigen/Core>

#include "gelex/logger.h"

namespace gelex
{

void constrain(Eigen::Ref<Eigen::VectorXd> varcmp, double y_variance)
{
    auto logger = logging::get();
    constexpr double constr_scale = 1e-6;
    const double limit = y_variance * constr_scale;

    Eigen::ArrayX<bool> mask = varcmp.array() < 0;
    auto num_constrained = mask.count();
    auto num_varcmp = varcmp.size();

    if (num_constrained == 0)
    {
        return;
    }
    if (num_constrained == varcmp.size())
    {
        logger->warn(
            "All variance components are constrained! The estimate is not "
            "reliable.");
        varcmp.fill(limit);
        return;
    }

    if (num_constrained > num_varcmp / 2)
    {
        logger->warn(
            "Half of the variance components are constrained! The estimate is "
            "not reliable.");
    }

    double delta = 0.0;
    for (int i = 0; i < num_varcmp; ++i)
    {
        if (mask[i])
        {
            delta += limit - varcmp[i];
            varcmp[i] = limit;
        }
    }

    delta /= static_cast<double>(num_varcmp - num_constrained);
    for (int i = 0; i < num_varcmp; ++i)
    {
        if (!mask[i] && varcmp[i] > delta)
        {
            varcmp[i] -= delta;
        }
    }
}
}  // namespace gelex
