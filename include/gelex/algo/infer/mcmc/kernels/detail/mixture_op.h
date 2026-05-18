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

#ifndef GELEX_ALGO_INFER_MCMC_KERNELS_DETAIL_MIXTURE_OP_H_
#define GELEX_ALGO_INFER_MCMC_KERNELS_DETAIL_MIXTURE_OP_H_

#include <variant>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/infra/stats/detail/var.h"

namespace gelex::mcmc::detail
{

constexpr int kMaxMixtureComponents = 5;

struct MixtureNormalPosteriors
{
    Eigen::Array<double, kMaxMixtureComponents, 1> means;
    Eigen::Array<double, kMaxMixtureComponents, 1> vars;
    Eigen::Array<double, kMaxMixtureComponents, 1> log_likelihoods;
};

template <typename StateT>
inline auto get_mixture(StateT& state) -> bayes::ComponentAllocation*
{
    if (!state.group)
    {
        return nullptr;
    }
    return std::get_if<bayes::ComponentAllocation>(&*state.group);
}

inline void compute_component_variances(
    bayes::ComponentAllocation& marker_assignment)
{
    for (Eigen::Index k = 0;
         k < static_cast<Eigen::Index>(marker_assignment.component_u.size());
         ++k)
    {
        marker_assignment.component_variance(k)
            = gelex::stats::detail::var(marker_assignment.component_u[k])(0);
    }
}

}  // namespace gelex::mcmc::detail

#endif  // GELEX_ALGO_INFER_MCMC_KERNELS_DETAIL_MIXTURE_OP_H_
