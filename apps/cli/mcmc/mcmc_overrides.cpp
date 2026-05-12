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

#include "mcmc_overrides.h"

#include <cmath>
#include <string_view>
#include <variant>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

namespace
{

auto check_side(
    const GeneticOverride& side,
    std::string_view pi_name,
    std::string_view mult_name) -> void
{
    if (side.proportions)
    {
        const auto& probs = *side.proportions;
        if ((probs.array() < 0.0).any())
        {
            throw GelexException(
                fmt::format("{} contains negative entries", pi_name));
        }
        if (std::abs(probs.sum() - 1.0) > 1e-9)
        {
            throw GelexException(
                fmt::format("{} must sum to 1, got {}", pi_name, probs.sum()));
        }
    }
    if (side.proportions && side.multiplier
        && side.proportions->size() != side.multiplier->size())
    {
        throw GelexException(
            fmt::format(
                "{} size ({}) must match {} size ({})",
                pi_name,
                side.proportions->size(),
                mult_name,
                side.multiplier->size()));
    }
}

auto apply_side(
    bayes::GeneticSpec& spec,
    std::optional<bayes::Mixture>& mixture,
    const GeneticOverride& side) -> void
{
    if (mixture && side.proportions)
    {
        mixture->proportions.init = *side.proportions;
    }
    if (mixture && side.multiplier)
    {
        if (auto* sm = std::get_if<bayes::ScaledMixture>(&mixture->strategy))
        {
            sm->multiplier = *side.multiplier;
        }
    }
    if (spec.mode == GeneticMode::D && side.positive_prob && spec.sign)
    {
        const double p = *side.positive_prob;
        spec.sign->init = Eigen::VectorXd{{p, 1.0 - p}};
    }
}

}  // namespace

auto validate_overrides(
    const bayes::LegacyBayesConfig& method,
    const MethodOverrides& overrides) -> void
{
    if (overrides.additive.positive_prob)
    {
        throw GelexException(
            "positive_prob is only valid for the dominance effect");
    }
    if (overrides.dominance.positive_prob)
    {
        const double p = *overrides.dominance.positive_prob;
        if (p < 0.0 || p > 1.0)
        {
            throw GelexException(
                fmt::format("positive_prob must be in [0, 1], got {}", p));
        }
        if (method.dominance != bayes::DominancePolicy::asymmetric)
        {
            throw GelexException(
                "positive_prob requires --asym dominance policy");
        }
    }
    check_side(overrides.additive, "pi", "multiplier");
    check_side(overrides.dominance, "dpi", "dmultiplier");
}

auto apply_overrides(
    bayes::LegacyBayesMethod& method,
    const MethodOverrides& overrides) -> void
{
    for (auto& prior : method.genetics)
    {
        if (auto* s = std::get_if<bayes::GeneticSpec>(&prior.spec))
        {
            apply_side(
                *s,
                prior.mixture,
                s->mode == GeneticMode::A ? overrides.additive
                                          : overrides.dominance);
        }
        else if (auto* j = std::get_if<bayes::JointSpec>(&prior.spec))
        {
            apply_side(j->additive, prior.mixture, overrides.additive);
            apply_side(j->dominance, prior.mixture, overrides.dominance);
        }
    }
}

}  // namespace gelex::cli
