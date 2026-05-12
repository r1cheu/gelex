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

#ifndef GELEX_CLI_MCMC_OVERRIDES_H_
#define GELEX_CLI_MCMC_OVERRIDES_H_

#include <optional>

#include <Eigen/Core>

#include "gelex/model/bayes/legacy_method.h"

namespace gelex::cli
{

struct GeneticOverride
{
    std::optional<Eigen::VectorXd> proportions;
    std::optional<Eigen::VectorXd> multiplier;
    std::optional<double> positive_prob;
};

struct MethodOverrides
{
    GeneticOverride additive;
    GeneticOverride dominance;
};

auto validate_overrides(
    const bayes::LegacyBayesConfig& method,
    const MethodOverrides& overrides) -> void;

auto apply_overrides(
    bayes::LegacyBayesMethod& method,
    const MethodOverrides& overrides) -> void;

}  // namespace gelex::cli

#endif  // GELEX_CLI_MCMC_OVERRIDES_H_
