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

#ifndef GELEX_MODEL_BAYES_PRIOR_NEW_H_
#define GELEX_MODEL_BAYES_PRIOR_NEW_H_

#include <cstdint>
#include <variant>

#include <Eigen/Core>

namespace gelex::bayes
{

struct ScaledInvChiSqPrior
{
    double nu{};
    double scale{};
};

struct DirichletPrior
{
    Eigen::VectorXi alpha;
};

enum class VarianceScope : std::uint8_t
{
    per_marker,
    per_block,
};

struct VarianceSpec
{
    VarianceScope scope{};
    double init{};
    ScaledInvChiSqPrior prior;
};

struct CategoricalSpec
{
    Eigen::VectorXd init;
    DirichletPrior prior;
    bool estimate = false;
};

struct SpikeSlab
{
};

struct ScaledMixture
{
    Eigen::VectorXd multiplier;
};

struct JointMixture
{
};

using VarianceStrategy = std::variant<SpikeSlab, ScaledMixture, JointMixture>;

struct Mixture
{
    VarianceStrategy strategy;
    CategoricalSpec proportions;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_PRIOR_NEW_H_
