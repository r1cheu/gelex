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

#include "gelex/model/bayes/bayes_policy.h"

#include <array>
#include <utility>

#include <Eigen/Core>

#include "gelex/model/bayes/prior.h"

namespace gelex::bayes
{

namespace
{

inline constexpr std::array kIndepShapes{
    AlgorithmShape::a_only,
    AlgorithmShape::d_only,
    AlgorithmShape::ad_independent,
};

const std::array<BayesPolicy, std::to_underlying(BayesBase::kCount)> kPolicies{{
    // A
    {
        .variance_scope = VarianceScope::per_marker,
        .mixture = std::nullopt,
        .supports_estimate_pi = false,
        .supports_asymmetric_dominance = false,
        .shapes = kIndepShapes,
    },
    // B
    {
        .variance_scope = VarianceScope::per_marker,
        .mixture = MixtureShape{
            Eigen::VectorXd{{0.99, 0.01}},
            SpikeSlab{},
        },
        .supports_estimate_pi = true,
        .supports_asymmetric_dominance = false,
        .shapes = kIndepShapes,
    },
    // C
    {
        .variance_scope = VarianceScope::per_block,
        .mixture = MixtureShape{
            Eigen::VectorXd{{0.99, 0.01}},
            SpikeSlab{},
        },
        .supports_estimate_pi = true,
        .supports_asymmetric_dominance = false,
        .shapes = kIndepShapes,
    },
    // R
    {
        .variance_scope = VarianceScope::per_block,
        .mixture = MixtureShape{
            Eigen::VectorXd{{0.99, 0.005, 0.001, 0.001, 0.001}},
            ScaledMixture{Eigen::VectorXd{{0.0, 0.001, 0.01, 0.1, 1.0}}},
        },
        .supports_estimate_pi = false,
        .supports_asymmetric_dominance = true,
        .shapes = kIndepShapes,
    },
    // RR
    {
        .variance_scope = VarianceScope::per_marker,
        .mixture = std::nullopt,
        .supports_estimate_pi = false,
        .supports_asymmetric_dominance = false,
        .shapes = kIndepShapes,
    },
}};

}  // namespace

auto policy_for(BayesBase base) -> const BayesPolicy&
{
    return kPolicies.at(std::to_underlying(base));
}

}  // namespace gelex::bayes
