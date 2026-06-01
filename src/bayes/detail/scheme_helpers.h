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

#ifndef GELEX_SRC_BAYES_DETAIL_SCHEME_HELPERS_H_
#define GELEX_SRC_BAYES_DETAIL_SCHEME_HELPERS_H_

#include <string_view>

#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/recipe_options.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{
class BayesModel;
}  // namespace gelex

namespace gelex::bayes::detail
{

auto effect_config(const BayesRecipeConfig& options, GeneticMode mode)
    -> const EffectConfig&;

auto heritability(GeneticMode mode, const EffectConfig& effect) -> double;

auto target_marker_variance(
    const BayesModel& model,
    GeneticMode mode,
    double h2,
    double active_marker_weight) -> double;

auto variance_parameter(double target) -> VarianceParameter;

auto mixture_proportion(const Simplex<double>& proportion, bool sampled)
    -> MixtureProportion;

auto scaled_active_marker_weight(
    const Simplex<double>& proportion,
    const ScaleMultiplier<double>& multiplier) -> double;

auto reject_joint_options(
    const BayesRecipeConfig& options,
    std::string_view scheme) -> void;

auto reject_proportion_options(
    const BayesRecipeConfig& options,
    std::string_view scheme) -> void;

auto reject_multiplier_options(
    const BayesRecipeConfig& options,
    std::string_view scheme) -> void;

auto reject_unpaired_proportion_multiplier(
    const BayesRecipeConfig& options,
    std::string_view scheme) -> void;

}  // namespace gelex::bayes::detail

#endif  // GELEX_SRC_BAYES_DETAIL_SCHEME_HELPERS_H_
