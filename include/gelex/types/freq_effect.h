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

#ifndef GELEX_TYPES_FREQ_EFFECT_H_
#define GELEX_TYPES_FREQ_EFFECT_H_

#include <cstdint>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <fmt/base.h>

#include "gelex/types/fixed_effects.h"

namespace gelex::freq
{

enum class GrmType : uint8_t
{
    A,
    D,
    AD,
    AA,
    DD,
    Unknown
};

struct RandomEffect
{
    std::string name;
    std::vector<std::string> levels;
    Eigen::MatrixXd K;
};

struct GeneticEffect
{
    GrmType type;
    Eigen::MatrixXd K;
};

struct FixedState
{
    explicit FixedState(const gelex::FixedEffect& effect);
    Eigen::VectorXd coeff;
    Eigen::VectorXd se;
};

struct RandomState
{
    explicit RandomState(const RandomEffect& effect);
    std::string name;
    Eigen::VectorXd blup;
    double variance{};
    double variance_se{};
};

struct GeneticState
{
    explicit GeneticState(const GeneticEffect& effect);
    GrmType type;
    Eigen::VectorXd ebv;
    double variance{};
    double variance_se{};
    double heritability{};
    double heritability_se{};
};

struct ResidualState
{
    double variance{};
    double variance_se{};
};

}  // namespace gelex::freq

namespace fmt
{
template <>
struct formatter<gelex::freq::GrmType> : formatter<string_view>
{
    auto format(gelex::freq::GrmType t, format_context& ctx) const
        -> format_context::iterator
    {
        string_view name = "unknown";
        switch (t)
        {
            case (gelex::freq::GrmType::A):
                name = "A";
                break;
            case (gelex::freq::GrmType::D):
                name = "D";
                break;
            case (gelex::freq::GrmType::AA):
                name = "AA";
                break;
            case (gelex::freq::GrmType::AD):
                name = "AD";
                break;
            case (gelex::freq::GrmType::DD):
                name = "DD";
                break;
            case (gelex::freq::GrmType::Unknown):
            default:
                break;
        }

        return formatter<string_view>::format(name, ctx);
    }
};
}  // namespace fmt

#endif  // GELEX_TYPES_FREQ_EFFECT_H_
