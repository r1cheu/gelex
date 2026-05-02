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

#ifndef GELEX_MODEL_FREQ_EFFECT_H_
#define GELEX_MODEL_FREQ_EFFECT_H_

#include <cstdint>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <fmt/base.h>

#include "gelex/types/fixed_effects.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::freq
{

struct RandomEffect
{
    std::string name;
    std::vector<std::string> levels;
    Eigen::MatrixXd K;
};

struct GeneticEffect
{
    GeneticMode type;
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
    GeneticMode type;
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

#endif  // GELEX_MODEL_FREQ_EFFECT_H_
