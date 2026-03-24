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

#ifndef GELEX_MODEL_BAYES_PRIOR_H_
#define GELEX_MODEL_BAYES_PRIOR_H_

#include <algorithm>
#include <optional>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class PriorSetConfig;
struct GeneticPriorConfig;

namespace bayes
{

struct GeneticEffect;

struct VariancePrior
{
    detail::ScaledInvChiSqParams param;
    double init;
    Eigen::Index size;
};

struct ProportionPrior
{
    Eigen::VectorXd init;
    bool estimate;
};

struct SignPrior
{
    double init_value;
};

// BayesA / BayesRR: all markers have non-zero effects
struct ContinuousPrior
{
    VariancePrior variance;
};

// BayesB / BayesC / BayesCpi: spike-and-slab (point mass at zero + slab)
struct SpikePrior
{
    VariancePrior variance;
    ProportionPrior proportion;
};

// BayesR: finite mixture of scaled normals
struct MixturePrior
{
    VariancePrior variance;
    ProportionPrior proportion;
    Eigen::VectorXd multiplier;
};

using MarkerPrior = std::variant<ContinuousPrior, SpikePrior, MixturePrior>;

struct GeneticPrior
{
    GeneticKind type;
    MarkerPrior marker;
    std::optional<SignPrior> sign;
};

using RandomPrior = VariancePrior;
using ResidualPrior = VariancePrior;

struct PriorSet
{
    static auto build(
        const PriorSetConfig& config,
        const std::vector<GeneticEffect>& genetics,
        std::size_t num_random_effects) -> PriorSet;

    [[nodiscard]] auto genetic(GeneticKind type) const -> const GeneticPrior*
    {
        auto it = std::ranges::find(genetics, type, &GeneticPrior::type);
        return it != genetics.end() ? &*it : nullptr;
    }

    [[nodiscard]] auto genetic(GeneticKind type) -> GeneticPrior*
    {
        auto it = std::ranges::find(genetics, type, &GeneticPrior::type);
        return it != genetics.end() ? &*it : nullptr;
    }

    std::vector<GeneticPrior> genetics;
    std::vector<RandomPrior> random;
    ResidualPrior residual;

   private:
    static auto build_genetic_prior(
        const GeneticEffect& effect,
        const GeneticPriorConfig& prior,
        const PriorSetConfig& config) -> GeneticPrior;
};

}  // namespace bayes
}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_PRIOR_H_
