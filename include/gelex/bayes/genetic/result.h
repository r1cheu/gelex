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
#include "gelex/bayes/genetic_family.h"

namespace gelex
{

namespace detail
{

template <VarianceLayout Kind>
using marker_variance_result_t = std::conditional_t<
    Kind == VarianceLayout::Pooled,
    ScalarPosteriorResult,
    EmptyPosteriorResult>;

template <MixtureWeightUpdate Update, typename Result>
using weight_result_t = std::conditional_t<
    Update == MixtureWeightUpdate::Enabled,
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

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
struct SpikeSlabPosteriorResult
{
    detail::marker_variance_result_t<Kind> variance;
    detail::weight_result_t<WeightUpdate, ScalarPosteriorResult> probability;
};

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
struct ScaledMixturePosteriorResult
{
    ScalarPosteriorResult variance;
    detail::weight_result_t<WeightUpdate, VectorPosteriorResult> probabilities;
    VectorPosteriorResult component_explained_variance;
};

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
struct JointSpikeSlabPosteriorResult
{
    detail::weight_result_t<WeightUpdate, VectorPosteriorResult> probabilities;
    VectorPosteriorResult component_explained_variance;
};

// Mirrors GeneticModeDraws: marker effects sit beside the family posterior.
template <typename FamilyResult>
struct GeneticModeResult
{
    VectorPosteriorResult coefficients;
    FamilyResult family_result;
};

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_RESULT_H_
