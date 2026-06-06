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

#include <array>
#include <vector>

#include "bayes/detail/scheme_helpers.h"
#include "gelex/bayes/genetic/half_normal_prior.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/recipe_options.h"
#include "gelex/bayes/scheme.h"
#include "gelex/exception.h"
#include "gelex/types/constrained_vector.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

BayesCDScheme::BayesCDScheme(const BayesRecipeOptions& options)
    : options_{options}
{
    if (options_.modes.size() != 2)
    {
        throw GelexException("CD requires --mode AD (both effects)");
    }
    if (options_.additive_proportion || options_.additive_proportion_update)
    {
        throw GelexException(
            "CD does not accept per-effect proportion override for A; "
            "use --jpi instead");
    }
    if (options_.dominance_proportion || options_.dominance_proportion_update)
    {
        throw GelexException(
            "CD does not accept per-effect proportion override for D; "
            "use --jpi instead");
    }
    if (options_.additive_multiplier)
    {
        throw GelexException(
            "CD does not accept per-effect multiplier override for A");
    }
    if (options_.dominance_multiplier)
    {
        throw GelexException(
            "CD does not accept per-effect multiplier override for D");
    }
    if (options_.joint_proportion && options_.joint_proportion->size() != 4)
    {
        throw GelexException(
            "CD requires joint_proportion to have exactly 4 components: "
            "both-off, A-only, D-only, both-on");
    }
}

auto BayesCDScheme::make_prior(const BayesModel& model) const
    -> std::vector<GeneticPrior>
{
    const Simplex<double> proportion = options_.joint_proportion.value_or(
        Simplex<double>{
            {0.99,
             0.003333333333333333,
             0.003333333333333333,
             0.003333333333333333}});
    // proportion layout: [both-off, A-only, D-only, both-on]
    const double active_a = proportion[1] + proportion[3];
    const double active_d = proportion[2] + proportion[3];
    const double dominance_positive_probability
        = options_.dominance_positive_probability
              ? options_.dominance_positive_probability->value()
              : 0.5;

    const double target_a = detail::target_marker_variance(
        model,
        GeneticMode::A,
        detail::heritability(options_, GeneticMode::A),
        active_a);
    const double target_d = detail::target_marker_variance(
        model,
        GeneticMode::D,
        detail::heritability(options_, GeneticMode::D),
        active_d);

    return std::vector<GeneticPrior>{
        JointGeneticPrior{JointHalfNormalMixturePrior{
            JointSharedMarkerVariance{std::array{
                SharedMarkerVariance{detail::variance_parameter(target_a)},
                SharedMarkerVariance{detail::variance_parameter(target_d)}}},
            detail::mixture_proportion(
                proportion, options_.joint_proportion_update.value_or(true)),
            ProbabilityParameter{
                dominance_positive_probability, BetaPrior{1.0, 1.0}}}}};
}

}  // namespace gelex::bayes
