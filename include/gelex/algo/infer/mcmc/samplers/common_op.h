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

#ifndef GELEX_ALGO_INFER_MCMC_SAMPLERS_COMMON_OP_H_
#define GELEX_ALGO_INFER_MCMC_SAMPLERS_COMMON_OP_H_

#include <cmath>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/infer/detail/marker_op.h"
#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/states.h"

namespace gelex::detail
{

constexpr int kMaxMixtureComponents = 10;

template <typename DerivedCol>
inline void update_component_u(
    std::vector<Eigen::VectorXd>& component_u,
    int old_index,
    double old_val,
    int new_index,
    double new_val,
    const Eigen::DenseBase<DerivedCol>& col)
{
    if (component_u.empty())
    {
        return;
    }

    if (old_index == new_index)
    {
        if (old_index > 0)
        {
            const double diff = new_val - old_val;
            if (std::abs(diff) > std::numeric_limits<double>::epsilon())
            {
                blas_daxpy(diff, col, component_u[old_index - 1]);
            }
        }
    }
    else
    {
        if (old_index > 0)
        {
            blas_daxpy(-old_val, col, component_u[old_index - 1]);
        }
        if (new_index > 0)
        {
            blas_daxpy(new_val, col, component_u[new_index - 1]);
        }
    }
}

inline void compute_component_variances(
    bayes::ComponentAllocation& marker_assignment)
{
    for (Eigen::Index k = 0;
         k < static_cast<Eigen::Index>(marker_assignment.component_u.size());
         ++k)
    {
        marker_assignment.component_variance(k)
            = detail::var(marker_assignment.component_u[k])(0);
    }
}

template <typename StateT>
inline void compute_component_variances(StateT& state)
{
    if (!state.group)
    {
        return;
    }
    auto* ma = std::get_if<bayes::ComponentAllocation>(&*state.group);
    if (ma && !ma->component_u.empty())
    {
        compute_component_variances(*ma);
    }
}

}  // namespace gelex::detail

#endif  // GELEX_ALGO_INFER_MCMC_SAMPLERS_COMMON_OP_H_
