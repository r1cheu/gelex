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

#ifndef GELEX_BAYES_LEGACY_PRIOR_H_
#define GELEX_BAYES_LEGACY_PRIOR_H_

#include <span>
#include <string_view>
#include <variant>
#include <vector>

#include "gelex/bayes/genetic/legacy_genetic_prior.h"
#include "gelex/bayes/parameter/values.h"

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
    auto visit(FieldVisitor& visitor) -> void
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
    auto visit(FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        parameter_.visit(visitor);
    }

   private:
    VarianceParameter parameter_;
};

using GeneticPrior = std::variant<SingleGeneticPrior, JointGeneticPrior>;

class LegacyBayesPrior
{
   public:
    static constexpr std::string_view name = "prior";

    LegacyBayesPrior(
        RandomPrior random,
        std::vector<GeneticPrior> genetics,
        ResidualPrior residual);

    LegacyBayesPrior(const LegacyBayesPrior&) = delete;
    LegacyBayesPrior(LegacyBayesPrior&&) noexcept = default;

    auto operator=(const LegacyBayesPrior&) -> LegacyBayesPrior& = delete;
    auto operator=(LegacyBayesPrior&&) noexcept -> LegacyBayesPrior& = default;

    ~LegacyBayesPrior() = default;

    auto random() const -> const RandomPrior& { return random_; }
    auto residual() const -> const ResidualPrior& { return residual_; }
    auto genetics() const -> std::span<const GeneticPrior> { return genetics_; }

    auto visit(FieldVisitor& visitor) -> void;

   private:
    static auto validate_genetics(const std::vector<GeneticPrior>& genetics)
        -> void;

    RandomPrior random_;
    std::vector<GeneticPrior> genetics_;
    ResidualPrior residual_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_LEGACY_PRIOR_H_
