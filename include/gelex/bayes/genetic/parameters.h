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

#include <array>
#include <ranges>
#include <string>
#include <string_view>

#include <Eigen/Core>

#include "gelex/bayes/parameter/values.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

class SharedMarkerVariance
{
   public:
    static constexpr std::string_view name = "shared_marker_variance";
    explicit SharedMarkerVariance(VarianceParameter parameter);

    auto parameter() -> VarianceParameter& { return parameter_; }
    auto parameter() const -> const VarianceParameter& { return parameter_; }
    auto visit(infra::FieldVisitor& visitor) -> void
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
    auto visit(infra::FieldVisitor& visitor) -> void
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
    auto visit(infra::FieldVisitor& visitor) -> void
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

class FixedProportion
{
   public:
    static constexpr std::string_view name = "fixed_mixture_proportion";
    explicit FixedProportion(Eigen::VectorXd value);

    auto value() const -> const Eigen::VectorXd& { return value_; }
    auto size() const -> Eigen::Index { return value_.size(); }
    auto visit(infra::FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        visitor.on("value", value_, FieldFlag::checkpoint | FieldFlag::report);
    }

   private:
    Eigen::VectorXd value_;
};

class SampledProportion
{
   public:
    static constexpr std::string_view name = "sampled_mixture_proportion";
    explicit SampledProportion(SimplexParameter parameter);

    auto parameter() const -> const SimplexParameter& { return parameter_; }
    auto size() const -> Eigen::Index { return parameter_.size(); }
    auto visit(infra::FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        parameter_.visit(visitor);
    }

   private:
    SimplexParameter parameter_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENETIC_PARAMETERS_H_
