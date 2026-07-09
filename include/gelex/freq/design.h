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

#include <Eigen/Core>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "gelex/types/fixed_designs.h"

namespace gelex::freq
{

enum class RandomKind : std::uint8_t
{
    Discrete,      // one-hot ZZ^T from a discrete factor
    Quantitative,  // linear kernel ZZ^T from a quantitative matrix
    Grm,           // genomic relationship matrix; chromosome-partitionable
};

struct RandomDesign
{
    std::string name;
    std::optional<std::vector<std::string>> levels;
    std::optional<Eigen::MatrixXd> Z;  // skip if identity
    Eigen::MatrixXd K;                 // Kernels
    RandomKind kind = RandomKind::Discrete;
};

struct RandomState
{
    explicit RandomState(const RandomDesign& design);
    RandomState() = default;
    Eigen::VectorXd blup;  // sample-level random predictions
    double variance{};
    double variance_se{};
    double variance_ratio{};
    double variance_ratio_se{};
};

struct FixedState
{
    explicit FixedState(const gelex::FixedDesign& design);
    Eigen::VectorXd coeffs;
    Eigen::VectorXd se;
};

struct ResidualState
{
    double variance{};
    double variance_se{};
};

}  // namespace gelex::freq

#endif  // GELEX_FREQ_DESIGN_H_
