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

#include <span>
#include <string_view>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/exception.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/model/bayes/state_capabilities.h"

namespace gelex::mcmc
{

inline auto make_variance_sampler(const bayes::VarianceParameter& parameter)
    -> stats::ScaledInvChi2Sampler<double>
{
    return stats::ScaledInvChi2Sampler<double>{
        parameter.prior().degrees_of_freedom(), parameter.prior().scale()};
}

inline auto require_marker_variances(
    const bayes::GeneticPrior& prior,
    std::string_view kernel_name) -> std::span<const bayes::MarkerVariance>
{
    const auto* cap = prior.query<bayes::MarkerVarianceCap>();
    if (cap == nullptr)
    {
        throw GelexException(
            fmt::format("{}: prior lacks variance capability", kernel_name));
    }
    return cap->variance();
}

inline auto require_variance_state(
    bayes::GeneticPriorState& state,
    std::string_view kernel_name) -> std::span<Eigen::VectorXd>
{
    auto* cap = state.query<bayes::VarianceStateCap>();
    if (cap == nullptr)
    {
        throw GelexException(
            fmt::format(
                "{}: prior state lacks variance capability", kernel_name));
    }
    return cap->variance();
}

inline auto require_multipliers(
    const bayes::GeneticPrior& prior,
    std::string_view kernel_name) -> std::span<const Eigen::VectorXd>
{
    const auto* cap = prior.query<bayes::MultiplierCap>();
    if (cap == nullptr)
    {
        throw GelexException(
            fmt::format("{}: prior lacks multiplier capability", kernel_name));
    }
    return cap->multiplier();
}

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_KERNELS_COMMON_H_
