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

#ifndef GELEX_BAYES_GENETIC_PARAMETERS_H_
#define GELEX_BAYES_GENETIC_PARAMETERS_H_

#include <Eigen/Core>
#include <array>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>

#include "gelex/bayes/parameter/values.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::bayes
{

class SharedMarkerVariance
{
   public:
    static constexpr std::string_view name = "shared_marker_variance";
    explicit SharedMarkerVariance(VarianceParameter parameter);

    auto parameter() -> VarianceParameter& { return parameter_; }
    auto parameter() const -> const VarianceParameter& { return parameter_; }
    auto visit(FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        parameter_.visit(visitor);
    }

   private:
    VarianceParameter parameter_;
};

class PerMarkerVariance
{
   public:
    static constexpr std::string_view name = "per_marker_variance";
    explicit PerMarkerVariance(VarianceParameter parameter);

    auto parameter() -> VarianceParameter& { return parameter_; }
    auto parameter() const -> const VarianceParameter& { return parameter_; }
    auto visit(FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        parameter_.visit(visitor);
    }

   private:
    VarianceParameter parameter_;
};

class JointSharedMarkerVariance
{
   public:
    static constexpr std::string_view name = "joint_shared_marker_variance";
    explicit JointSharedMarkerVariance(
        std::array<SharedMarkerVariance, 2> variances);

    auto variance(GeneticMode mode) -> SharedMarkerVariance&
    {
        return variances_[std::to_underlying(mode)];
    }
    auto variance(GeneticMode mode) const -> const SharedMarkerVariance&
    {
        return variances_[std::to_underlying(mode)];
    }
    auto visit(FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        constexpr std::array modes{GeneticMode::A, GeneticMode::D};
        for (auto [i, mode_value] : std::views::enumerate(modes))
        {
            auto mode_scope = visitor.scope(std::to_string(i));
            auto mode = mode_value;
            visitor.on("mode", mode, FieldFlag::checkpoint | FieldFlag::report);
            variances_[i].visit(visitor);
        }
    }

   private:
    std::array<SharedMarkerVariance, 2> variances_;
};

class MixtureProportion
{
   public:
    static constexpr std::string_view name = "mixture_proportion";
    explicit MixtureProportion(Eigen::VectorXd initial_value);
    explicit MixtureProportion(SimplexParameter parameter);

    auto initial_value() const -> const Eigen::VectorXd&
    {
        return initial_value_;
    }
    auto is_sampled() const -> bool { return prior_.has_value(); }
    auto prior() const -> const std::optional<DirichletPrior>&
    {
        return prior_;
    }
    auto size() const -> Eigen::Index { return initial_value_.size(); }
    auto visit(FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        visitor.on(
            "initial_value",
            initial_value_,
            FieldFlag::checkpoint | FieldFlag::report);
        if (prior_)
        {
            prior_->visit(visitor);
        }
    }

   private:
    Eigen::VectorXd initial_value_;
    std::optional<DirichletPrior> prior_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENETIC_PARAMETERS_H_
