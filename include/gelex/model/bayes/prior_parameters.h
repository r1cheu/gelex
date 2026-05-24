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

#ifndef GELEX_MODEL_BAYES_PRIOR_PARAMETERS_H_
#define GELEX_MODEL_BAYES_PRIOR_PARAMETERS_H_

#include <cstdint>
#include <string_view>
#include <utility>

#include <Eigen/Core>

#include "gelex/infra/field_visitor.h"

namespace gelex::bayes
{

class ScaledInvChiSqPrior
{
   public:
    static constexpr std::string_view name = "scaled_inv_chi2";
    ScaledInvChiSqPrior() = default;
    ScaledInvChiSqPrior(double degrees_of_freedom, double scale);

    auto degrees_of_freedom() const -> double { return degrees_of_freedom_; }
    auto scale() const -> double { return scale_; }
    auto visit(infra::FieldVisitor& visitor) const -> void
    {
        auto scope = visitor.scope(name);
        visitor.on(
            "degrees_of_freedom", degrees_of_freedom_, FieldFlag::source);
        visitor.on("scale", scale_, FieldFlag::source);
    }

   private:
    double degrees_of_freedom_{-2};
    double scale_{0};
};

class DirichletPrior
{
   public:
    static constexpr std::string_view name = "dirichlet";
    explicit DirichletPrior(Eigen::VectorXd concentration);

    auto concentration() const -> const Eigen::VectorXd&
    {
        return concentration_;
    }
    auto size() const -> Eigen::Index { return concentration_.size(); }
    auto visit(infra::FieldVisitor& visitor) const -> void
    {
        auto scope = visitor.scope(name);
        visitor.on("concentration", concentration_, FieldFlag::source);
    }

   private:
    Eigen::VectorXd concentration_;
};

class VarianceParameter
{
   public:
    static constexpr std::string_view name = "variance";
    VarianceParameter(double initial_value, ScaledInvChiSqPrior prior);

    auto initial_value() const -> double { return initial_value_; }
    auto prior() const -> const ScaledInvChiSqPrior& { return prior_; }
    auto visit(infra::FieldVisitor& visitor) const -> void
    {
        auto scope = visitor.scope(name);
        visitor.on("initial_value", initial_value_, FieldFlag::source);
        prior_.visit(visitor);
    }

   private:
    double initial_value_;
    ScaledInvChiSqPrior prior_;
};

class SimplexParameter
{
   public:
    static constexpr std::string_view name = "simplex";
    SimplexParameter(Eigen::VectorXd initial_value, DirichletPrior prior);

    auto initial_value() const -> const Eigen::VectorXd&
    {
        return initial_value_;
    }
    auto prior() const -> const DirichletPrior& { return prior_; }
    auto size() const -> Eigen::Index { return initial_value_.size(); }
    auto visit(infra::FieldVisitor& visitor) const -> void
    {
        auto scope = visitor.scope(name);
        visitor.on("initial_value", initial_value_, FieldFlag::source);
        prior_.visit(visitor);
    }

   private:
    Eigen::VectorXd initial_value_;
    DirichletPrior prior_;
};

enum class MarkerVarianceLayout : std::uint8_t
{
    shared,
    per_marker,
};

class MarkerVariance
{
   public:
    static constexpr std::string_view name = "marker_variance";
    MarkerVariance(MarkerVarianceLayout layout, VarianceParameter parameter);

    auto layout() const -> MarkerVarianceLayout { return layout_; }
    auto parameter() const -> const VarianceParameter& { return parameter_; }
    auto marker_variance_size(Eigen::Index num_markers) const -> Eigen::Index
    {
        switch (layout_)
        {
            case MarkerVarianceLayout::shared:
                return 1;
            case MarkerVarianceLayout::per_marker:
                return num_markers;
        }
        std::unreachable();
    }
    auto visit(infra::FieldVisitor& visitor) const -> void
    {
        auto scope = visitor.scope(name);
        visitor.on("layout", layout_, FieldFlag::source | FieldFlag::metadata);
        parameter_.visit(visitor);
    }

   private:
    MarkerVarianceLayout layout_{MarkerVarianceLayout::shared};
    VarianceParameter parameter_;
};

enum class UpdatePolicy : std::uint8_t
{
    fixed,
    sampled,
};

class MixtureProportion
{
   public:
    static constexpr std::string_view name = "mixture_proportion";
    MixtureProportion(SimplexParameter parameter, UpdatePolicy update);

    auto parameter() const -> const SimplexParameter& { return parameter_; }
    auto update() const -> UpdatePolicy { return update_; }

    auto size() const -> Eigen::Index { return parameter_.size(); }
    auto sampled() const -> bool { return update_ == UpdatePolicy::sampled; }
    auto visit(infra::FieldVisitor& visitor) const -> void
    {
        auto scope = visitor.scope(name);
        parameter_.visit(visitor);
        visitor.on("update", update_, FieldFlag::source | FieldFlag::metadata);
    }

   private:
    SimplexParameter parameter_;
    UpdatePolicy update_{UpdatePolicy::fixed};
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_PRIOR_PARAMETERS_H_
