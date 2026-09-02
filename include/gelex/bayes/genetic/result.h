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

#ifndef GELEX_BAYES_GENETIC_RESULT_H_
#define GELEX_BAYES_GENETIC_RESULT_H_

#include <cstddef>
#include <type_traits>

#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

namespace detail
{

template <VarianceLayout Kind>
using marker_variance_result_t = std::conditional_t<
    Kind == VarianceLayout::Pooled,
    ScalarPosteriorResult,
    EmptyPosteriorResult>;

template <UpdatePolicy Policy, typename Result>
using policy_result_t = std::conditional_t<
    Policy == UpdatePolicy::Sampled,
    Result,
    EmptyPosteriorResult>;

}  // namespace detail

template <VarianceLayout Kind>
struct GaussianPosteriorResult
{
    detail::marker_variance_result_t<Kind> variance;
};

template <HalfNormalAsymmetry Axis>
struct HalfNormalPosteriorResult;

template <>
struct HalfNormalPosteriorResult<HalfNormalAsymmetry::Count>
{
    ScalarPosteriorResult variance;
    ScalarPosteriorResult positive_probability;
};

template <>
struct HalfNormalPosteriorResult<HalfNormalAsymmetry::Magnitude>
{
    VectorPosteriorResult variances;
};

template <VarianceLayout Kind, UpdatePolicy ProbabilityUpdate>
struct SpikeSlabPosteriorResult
{
    detail::marker_variance_result_t<Kind> variance;
    detail::policy_result_t<ProbabilityUpdate, ScalarPosteriorResult>
        probability;
};

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
struct ScaledMixturePosteriorResult
{
    ScalarPosteriorResult variance;
    detail::policy_result_t<ProbabilitiesUpdate, VectorPosteriorResult>
        probabilities;
    VectorPosteriorResult component_explained_variance;
};

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
struct JointSpikeSlabPosteriorResult
{
    detail::policy_result_t<ProbabilitiesUpdate, VectorPosteriorResult>
        probabilities;
    VectorPosteriorResult component_explained_variance;
};

// Mirrors GeneticModeDraws: marker effects sit beside the family posterior.
template <typename FamilyResult>
struct GeneticModeResult
{
    VectorPosteriorResult coefficients;
    FamilyResult family_result;
};

namespace detail
{

template <typename Prior>
struct FamilyResultType;

template <VarianceLayout Kind>
struct FamilyResultType<GaussianPrior<Kind>>
{
    using type = GaussianPosteriorResult<Kind>;
};

template <HalfNormalAsymmetry Axis>
struct FamilyResultType<HalfNormalPrior<Axis>>
{
    using type = HalfNormalPosteriorResult<Axis>;
};

template <VarianceLayout Kind, UpdatePolicy ProbabilityUpdate>
struct FamilyResultType<SpikeSlabPrior<Kind, ProbabilityUpdate>>
{
    using type = SpikeSlabPosteriorResult<Kind, ProbabilityUpdate>;
};

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
struct FamilyResultType<ScaledMixturePrior<ClassCount, ProbabilitiesUpdate>>
{
    using type = ScaledMixturePosteriorResult<ClassCount, ProbabilitiesUpdate>;
};

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
struct FamilyResultType<JointSpikeSlabPrior<ClassCount, ProbabilitiesUpdate>>
{
    using type = JointSpikeSlabPosteriorResult<ClassCount, ProbabilitiesUpdate>;
};

template <typename Prior>
using family_result_t = typename FamilyResultType<Prior>::type;

template <typename GeneticPrior>
struct GeneticResultType;

template <GeneticModeSet Modes, typename... Priors>
struct GeneticResultType<ModeValues<Modes, Priors...>>
{
    using type
        = ModeValues<Modes, GeneticModeResult<family_result_t<Priors>>...>;
};

template <typename ModeValuesType, typename JointPrior>
struct GeneticResultType<JointModeValues<ModeValuesType, JointPrior>>
{
    using type = JointModeValues<
        typename GeneticResultType<ModeValuesType>::type,
        family_result_t<JointPrior>>;
};

template <typename GeneticPrior>
using genetic_result_t = typename GeneticResultType<GeneticPrior>::type;

}  // namespace detail

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_RESULT_H_
