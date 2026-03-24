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

#ifndef GELEX_MODEL_BAYES_STATES_H_
#define GELEX_MODEL_BAYES_STATES_H_

#include <cstdint>
#include <optional>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/types/fixed_effects.h"

namespace gelex::bayes
{

using TrackerVector = Eigen::VectorX<int8_t>;

struct Assignment
{
    Assignment(Eigen::Index num_markers, const Eigen::VectorXd& init_proportion)
        : tracker(TrackerVector::Zero(num_markers)),
          proportion(init_proportion),
          count(Eigen::VectorXi::Zero(init_proportion.size()))
    {
    }

    TrackerVector tracker;
    Eigen::VectorXd proportion;
    Eigen::VectorXi count;
};

struct ComponentAllocation
{
    ComponentAllocation(
        Eigen::Index num_markers,
        Eigen::Index num_samples,
        const Eigen::VectorXd& init_proportion)
        : assignment(num_markers, init_proportion),
          component_u(init_proportion.size() - 1),
          component_variance(Eigen::VectorXd::Zero(init_proportion.size() - 1))
    {
        for (auto& vec : component_u)
        {
            vec = Eigen::VectorXd::Zero(num_samples);
        }
    }

    Assignment assignment;
    std::vector<Eigen::VectorXd> component_u;
    Eigen::VectorXd component_variance;
};

using MarkerAllocation = std::variant<Assignment, ComponentAllocation>;

struct FixedState
{
    explicit FixedState(const FixedEffect& effect)
        : coeffs(Eigen::VectorXd::Zero(effect.X.cols())) {};
    Eigen::VectorXd coeffs;
};

struct RandomState
{
    RandomState(const RandomEffect& effect, const RandomPrior& prior)
        : coeffs(Eigen::VectorXd::Zero(effect.X.cols())), variance{prior.init}
    {
    }

    Eigen::VectorXd coeffs;
    double variance{0.0};
};

struct GeneticState
{
    GeneticState(const GeneticEffect& effect, const GeneticPrior& prior)
        : type(effect.type),
          coeffs(Eigen::VectorXd::Zero(bayes::get_cols(effect.X))),
          u(Eigen::VectorXd::Zero(bayes::get_rows(effect.X)))
    {
        auto num_markers = bayes::get_cols(effect.X);
        auto num_samples = bayes::get_rows(effect.X);

        std::visit(
            [&](const auto& p)
            {
                using T = std::decay_t<decltype(p)>;
                marker_variance = Eigen::VectorXd::Constant(
                    p.variance.size, p.variance.init);

                if constexpr (std::is_same_v<T, SpikePrior>)
                {
                    group.emplace(Assignment(num_markers, p.proportion.init));
                }
                else if constexpr (std::is_same_v<T, MixturePrior>)
                {
                    group.emplace(ComponentAllocation(
                        num_markers, num_samples, p.proportion.init));
                }
            },
            prior.marker);

        if (prior.sign)
        {
            Eigen::Vector3d sign_proportion{
                {0.0, prior.sign->init_value, 1.0 - prior.sign->init_value}};
            sign.emplace(num_markers, sign_proportion);
        }
    }

    GeneticKind type;
    Eigen::VectorXd coeffs;
    Eigen::VectorXd u;

    double variance{};
    double heritability{};
    Eigen::VectorXd marker_variance;

    std::optional<MarkerAllocation> group;
    std::optional<Assignment> sign;
};

struct ResidualState
{
    Eigen::VectorXd y_adj;
    double variance{0.0};
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_STATES_H_
