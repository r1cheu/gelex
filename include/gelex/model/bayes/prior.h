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
    stats::detail::ScaledInvChiSqParams param;
    double init{};
    Eigen::Index size{};
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
    GeneticMode type;
    MarkerPrior marker;
    std::optional<SignPrior> sign;
};

using RandomPrior = VariancePrior;
using ResidualPrior = VariancePrior;

class Priors
{
   public:
    Priors(
        const PriorSetConfig& config,
        const std::vector<GeneticEffect>& genetics,
        std::size_t num_random_effects);

    Priors(
        std::vector<GeneticPrior> genetics,
        std::vector<RandomPrior> random,
        ResidualPrior residual);

    [[nodiscard]] auto genetic(GeneticMode type) const -> const GeneticPrior*
    {
        auto it = std::ranges::find(genetics_, type, &GeneticPrior::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    [[nodiscard]] auto genetics() const -> const std::vector<GeneticPrior>&
    {
        return genetics_;
    }

    [[nodiscard]] auto random() const -> const std::vector<RandomPrior>&
    {
        return random_;
    }

    [[nodiscard]] auto residual() const -> const ResidualPrior&
    {
        return residual_;
    }

   private:
    std::vector<GeneticPrior> genetics_;
    std::vector<RandomPrior> random_;
    ResidualPrior residual_;

    static auto build_genetic_prior(
        const GeneticEffect& effect,
        const GeneticPriorConfig& prior,
        const PriorSetConfig& config) -> GeneticPrior;
};

}  // namespace bayes
}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_PRIOR_H_
