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

#ifndef GELEX_ALGO_INFER_MCMC_KERNELS_COMMON_H_
#define GELEX_ALGO_INFER_MCMC_KERNELS_COMMON_H_

#include <string_view>
#include <variant>

#include <fmt/format.h>

#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/exception.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/prior.h"

namespace gelex::mcmc
{

inline auto make_variance_sampler(const bayes::OldVarianceSpec& spec)
    -> stats::ScaledInvChi2Sampler<double>
{
    return stats::ScaledInvChi2Sampler<double>{
        spec.prior.degrees_of_freedom(), spec.prior.scale()};
}

template <typename Allocation>
auto unpack_marker_allocation(
    bayes::GeneticState& state,
    std::string_view kernel_name) -> Allocation&
{
    if (!state.group.has_value())
    {
        throw GelexException(
            fmt::format(
                "{}: state.group is empty (expected MarkerAllocation)",
                kernel_name));
    }
    if (auto* p = std::get_if<Allocation>(&*state.group); p != nullptr)
    {
        return *p;
    }
    throw GelexException(
        fmt::format(
            "{}: marker allocation variant mismatch (actual_index={})",
            kernel_name,
            state.group->index()));
}

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_KERNELS_COMMON_H_
