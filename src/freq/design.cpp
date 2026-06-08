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

#include "gelex/freq/design.h"

#include <span>
#include <string_view>

#include <Eigen/Core>

#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/types/fixed_designs.h"

namespace gelex::freq
{

FixedState::FixedState(const gelex::FixedDesign& design)
    : coeffs(Eigen::VectorXd::Zero(design.X.cols())),
      se(Eigen::VectorXd::Zero(design.X.cols()))
{
}

auto FixedState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("coeffs", coeffs, FieldFlag::report);
    visitor.on("se", se, FieldFlag::report);
}

RandomState::RandomState(const RandomDesign& design)
    : blup(Eigen::VectorXd::Zero(design.K.rows()))
{
}

auto RandomDesign::visit(infra::FieldVisitor& visitor) const -> void
{
    auto scope = visitor.scope(name);
    visitor.on("name", term_name, FieldFlag::report);
    if (levels)
    {
        visitor.on(
            "levels",
            std::span<const std::string>{levels->data(), levels->size()},
            FieldFlag::report);
    }
}

auto RandomState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("var", variance, FieldFlag::report);
    visitor.on("var_se", variance_se, FieldFlag::report);
    visitor.on("var_ratio", variance_ratio, FieldFlag::report);
    visitor.on("var_ratio_se", variance_ratio_se, FieldFlag::report);
}

auto ResidualState::visit(infra::FieldVisitor& visitor) -> void
{
    auto scope = visitor.scope(name);
    visitor.on("var", variance, FieldFlag::report);
    visitor.on("var_se", variance_se, FieldFlag::report);
}

}  // namespace gelex::freq
