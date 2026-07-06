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

#ifndef GELEX_BAYES_PARAMETER_DISTRIBUTIONS_H_
#define GELEX_BAYES_PARAMETER_DISTRIBUTIONS_H_

#include <string_view>

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
    auto visit(FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        visitor.on(
            "degrees_of_freedom",
            degrees_of_freedom_,
            FieldFlag::checkpoint | FieldFlag::report);
        visitor.on("scale", scale_, FieldFlag::checkpoint | FieldFlag::report);
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
    auto visit(FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        visitor.on(
            "concentration",
            concentration_,
            FieldFlag::checkpoint | FieldFlag::report);
    }

   private:
    Eigen::VectorXd concentration_;
};

class BetaPrior
{
   public:
    static constexpr std::string_view name = "beta";
    BetaPrior(double alpha, double beta);

    auto alpha() const -> double { return alpha_; }
    auto beta() const -> double { return beta_; }
    auto visit(FieldVisitor& visitor) -> void
    {
        auto scope = visitor.scope(name);
        visitor.on("alpha", alpha_, FieldFlag::checkpoint | FieldFlag::report);
        visitor.on("beta", beta_, FieldFlag::checkpoint | FieldFlag::report);
    }

   private:
    double alpha_;
    double beta_;
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_PARAMETER_DISTRIBUTIONS_H_
