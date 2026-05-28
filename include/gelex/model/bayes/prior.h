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
#include <span>
#include <string_view>
#include <variant>
#include <vector>

#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior_parameters.h"

namespace gelex::bayes
{

class RandomPrior
{
   public:
    static constexpr std::string_view name = "random";
    RandomPrior(double initial_value, ScaledInvChiSqPrior prior)
        : parameter_(initial_value, prior)
    {
    }

    explicit RandomPrior(VarianceParameter parameter) : parameter_(parameter) {}

    auto initial_value() const -> double { return parameter_.initial_value(); }
    auto prior() const -> const ScaledInvChiSqPrior&
    {
        return parameter_.prior();
    }
    auto visit(infra::FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        parameter_.visit(visitor);
    }

   private:
    VarianceParameter parameter_;
};

class ResidualPrior
{
   public:
    static constexpr std::string_view name = "residual";
    ResidualPrior(double initial_value, ScaledInvChiSqPrior prior)
        : parameter_(initial_value, prior)
    {
    }

    explicit ResidualPrior(VarianceParameter parameter) : parameter_(parameter)
    {
    }

    auto initial_value() const -> double { return parameter_.initial_value(); }
    auto prior() const -> const ScaledInvChiSqPrior&
    {
        return parameter_.prior();
    }
    auto visit(infra::FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        parameter_.visit(visitor);
    }

   private:
    VarianceParameter parameter_;
};

using GeneticPriorBlock = std::variant<
    std::unique_ptr<SingleGeneticPrior>,
    std::unique_ptr<JointGeneticPrior>>;

class BayesPrior
{
   public:
    static constexpr std::string_view name = "prior";

    BayesPrior(
        RandomPrior random,
        std::vector<GeneticPriorBlock> genetics,
        ResidualPrior residual);

    BayesPrior(const BayesPrior&) = delete;
    BayesPrior(BayesPrior&&) noexcept = default;

    auto operator=(const BayesPrior&) -> BayesPrior& = delete;
    auto operator=(BayesPrior&&) noexcept -> BayesPrior& = default;

    ~BayesPrior() = default;

    auto random() const -> const RandomPrior& { return random_; }
    auto residual() const -> const ResidualPrior& { return residual_; }
    auto genetics() const -> std::span<const GeneticPriorBlock>
    {
        return genetics_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void;

   private:
    static auto validate_genetics(
        const std::vector<GeneticPriorBlock>& genetics) -> void;

    RandomPrior random_;
    std::vector<GeneticPriorBlock> genetics_;
    ResidualPrior residual_;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_PRIOR_H_
