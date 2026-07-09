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

#ifndef GELEX_BAYES_PARAMETER_VALUES_H_
#define GELEX_BAYES_PARAMETER_VALUES_H_

#include <Eigen/Core>
#include <string_view>

#include "gelex/bayes/parameter/distributions.h"
#include "gelex/infra/field_visitor.h"

namespace gelex::bayes
{

class VarianceParameter
{
   public:
    static constexpr std::string_view name = "variance_param";
    VarianceParameter(double initial_value, ScaledInvChiSqPrior prior);

    auto initial_value() const -> double { return initial_value_; }
    auto prior() const -> const ScaledInvChiSqPrior& { return prior_; }
    auto visit(FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        visitor.on(
            "initial_value",
            initial_value_,
            FieldFlag::checkpoint | FieldFlag::report);
        prior_.visit(visitor);
    }

   private:
    double initial_value_;
    ScaledInvChiSqPrior prior_;
};

class SimplexParameter
{
   public:
    static constexpr std::string_view name = "simplex_param";
    SimplexParameter(Eigen::VectorXd initial_value, DirichletPrior prior);

    auto initial_value() const -> const Eigen::VectorXd&
    {
        return initial_value_;
    }
    auto prior() const -> const DirichletPrior& { return prior_; }
    auto size() const -> Eigen::Index { return initial_value_.size(); }
    auto visit(FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        visitor.on(
            "initial_value",
            initial_value_,
            FieldFlag::checkpoint | FieldFlag::report);
        prior_.visit(visitor);
    }

   private:
    Eigen::VectorXd initial_value_;
    DirichletPrior prior_;
};

class ProbabilityParameter
{
   public:
    static constexpr std::string_view name = "probability_param";
    ProbabilityParameter(double initial_value, BetaPrior prior);

    auto initial_value() const -> double { return initial_value_; }
    auto prior() const -> const BetaPrior& { return prior_; }
    auto visit(FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        visitor.on(
            "initial_value",
            initial_value_,
            FieldFlag::checkpoint | FieldFlag::report);
        prior_.visit(visitor);
    }

   private:
    double initial_value_;
    BetaPrior prior_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_PARAMETER_VALUES_H_
