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

#ifndef GELEX_FREQ_DESIGN_H_
#define GELEX_FREQ_DESIGN_H_

#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include <Eigen/Core>

#include <fmt/base.h>

#include "gelex/types/fixed_designs.h"

namespace gelex::infra
{
class FieldVisitor;
}

namespace gelex::freq
{

struct RandomDesign
{
    static constexpr std::string_view name{"random"};
    std::string term_name;
    std::optional<std::vector<std::string>> levels;
    std::optional<Eigen::MatrixXd> Z;  // skip if identity
    Eigen::MatrixXd K;                 // Kernels

    auto visit(infra::FieldVisitor& visitor) const -> void;
};

struct RandomState
{
    static constexpr std::string_view name{"random"};
    explicit RandomState(const RandomDesign& design);
    RandomState() = default;
    Eigen::VectorXd blup;  // sample-level random predictions
    double variance{};
    double variance_se{};
    double variance_ratio{};
    double variance_ratio_se{};

    auto visit(infra::FieldVisitor& visitor) -> void;
};

struct FixedState
{
    static constexpr std::string_view name{"fixed"};

    explicit FixedState(const gelex::FixedDesign& design);
    Eigen::VectorXd coeffs;
    Eigen::VectorXd se;

    auto visit(infra::FieldVisitor& visitor) -> void;
};

struct ResidualState
{
    static constexpr std::string_view name{"residual"};
    double variance{};
    double variance_se{};

    auto visit(infra::FieldVisitor& visitor) -> void;
};

}  // namespace gelex::freq

#endif  // GELEX_FREQ_DESIGN_H_
