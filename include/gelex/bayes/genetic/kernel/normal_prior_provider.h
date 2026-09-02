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

#ifndef GELEX_BAYES_GENETIC_KERNEL_NORMAL_PRIOR_PROVIDER_H_
#define GELEX_BAYES_GENETIC_KERNEL_NORMAL_PRIOR_PROVIDER_H_

#include <Eigen/Core>

#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/stats/quadratic_log_kernel.h"

namespace gelex::detail
{

template <VarianceLayout Kind>
[[nodiscard]] auto make_normal_prior_provider(
    const marker_variance_state_t<Kind>& variance)
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return [prior = make_normal_prior(variance)](
                   Eigen::Index) -> const QuadraticLogKernel& { return prior; };
    }
    else
    {
        return [variance = &variance](Eigen::Index marker)
        { return make_normal_prior((*variance)(marker)); };
    }
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_KERNEL_NORMAL_PRIOR_PROVIDER_H_
