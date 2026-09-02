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

#ifndef GELEX_BAYES_GENETIC_LEGACY_GENETIC_PRIOR_H_
#define GELEX_BAYES_GENETIC_LEGACY_GENETIC_PRIOR_H_

#include <Eigen/Core>
#include <string_view>
#include <variant>

#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/half_normal_prior.h"
#include "gelex/bayes/genetic/prior_state.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{
class FieldVisitor;
}

namespace gelex::bayes
{

inline constexpr std::string_view single_genetic_prior_name = "single";
inline constexpr std::string_view joint_genetic_prior_name = "joint";

using SingleGeneticPrior = std::variant<
    SingleSharedGaussianPrior,
    SinglePerMarkerGaussianPrior,
    SingleSharedSpikeSlabGaussianPrior,
    SinglePerMarkerSpikeSlabGaussianPrior,
    SingleScaledMixtureGaussianPrior>;

using JointGeneticPrior
    = std::variant<JointGaussianMixturePrior, JointHalfNormalMixturePrior>;

auto mode(const SingleGeneticPrior& prior) -> GeneticMode;

auto visit(SingleGeneticPrior& prior, FieldVisitor& visitor) -> void;
auto visit(JointGeneticPrior& prior, FieldVisitor& visitor) -> void;

auto make_state(
    const SingleGeneticPrior& prior,
    Eigen::Index num_markers,
    Eigen::Index num_individuals) -> SingleGeneticPriorState;
auto make_state(
    const JointGeneticPrior& prior,
    Eigen::Index num_markers,
    Eigen::Index num_individuals) -> JointGeneticPriorState;

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENETIC_LEGACY_GENETIC_PRIOR_H_
