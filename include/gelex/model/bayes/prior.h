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

#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior_specs.h"

namespace gelex::bayes
{

struct OldDirichletPrior
{
    Eigen::VectorXi concentration;
};

struct OldVarianceSpec
{
    MarkerVarianceScope scope{};
    double init{};
    ScaledInvChiSqPrior prior;

    static auto make(double phenotype_variance) -> OldVarianceSpec;
    static auto make(double init, MarkerVarianceScope scope) -> OldVarianceSpec;
};

struct CategoricalSpec
{
    Eigen::VectorXd init;
    OldDirichletPrior prior;
    bool estimate = false;
};

struct SpikeSlab
{
};

struct ScaledMixture
{
    Eigen::VectorXd multiplier;
};

struct JointMixture
{
};

using VarianceStrategy = std::variant<SpikeSlab, ScaledMixture, JointMixture>;

struct BayesPolicy;

struct Mixture
{
    VarianceStrategy strategy;
    CategoricalSpec proportions;

    static auto make(const BayesPolicy&, bool estimate_pi)
        -> std::optional<Mixture>;
};

struct RandomEffectPrior
{
    std::string name;
    VarianceSpec variance;
};

class BayesPrior
{
   public:
    BayesPrior(
        std::vector<RandomEffectPrior> randoms,
        std::vector<std::unique_ptr<GeneticPrior>> genetics,
        VarianceSpec residual);

    BayesPrior(const BayesPrior&) = delete;
    BayesPrior(BayesPrior&&) noexcept = default;

    auto operator=(const BayesPrior&) -> BayesPrior& = delete;
    auto operator=(BayesPrior&&) noexcept -> BayesPrior& = default;

    ~BayesPrior() = default;

    auto randoms() const -> std::span<const RandomEffectPrior>
    {
        return randoms_;
    }
    auto residual() const -> const VarianceSpec& { return residual_; }

    auto genetics() const -> decltype(auto)
    {
        return genetics_
               | std::views::transform(
                   [](const auto& prior) -> const GeneticPrior&
                   { return *prior; });
    }

    auto random(std::string_view name) const -> const RandomEffectPrior*;

   private:
    std::vector<RandomEffectPrior> randoms_;
    std::vector<std::unique_ptr<GeneticPrior>> genetics_;
    VarianceSpec residual_;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_PRIOR_H_
