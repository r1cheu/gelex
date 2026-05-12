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

#ifndef GELEX_ALGO_INFER_MCMC_STATE_H_
#define GELEX_ALGO_INFER_MCMC_STATE_H_

#include <cstdint>
#include <optional>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/types/fixed_effects.h"
#include "gelex/types/genetic_effect_type.h"

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
        count(0) = static_cast<int>(num_markers);
    }

    Assignment(
        TrackerVector tracker,
        Eigen::VectorXd proportion,
        Eigen::VectorXi count)
        : tracker(std::move(tracker)),
          proportion(std::move(proportion)),
          count(std::move(count))
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

    ComponentAllocation(
        Assignment assignment,
        std::vector<Eigen::VectorXd> component_u,
        Eigen::VectorXd component_variance)
        : assignment(std::move(assignment)),
          component_u(std::move(component_u)),
          component_variance(std::move(component_variance))
    {
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
    explicit FixedState(Eigen::VectorXd coeffs) : coeffs(std::move(coeffs)) {}
    Eigen::VectorXd coeffs;
};

struct RandomState
{
    RandomState(const RandomEffect& effect, const OldVarianceSpec& spec)
        : coeffs(Eigen::VectorXd::Zero(effect.X.cols())), variance{spec.init}
    {
    }

    RandomState(Eigen::VectorXd coeffs, double variance)
        : coeffs(std::move(coeffs)), variance{variance}
    {
    }

    Eigen::VectorXd coeffs;
    double variance{0.0};
};

struct GeneticState
{
    GeneticState(
        GeneticMode type,
        Eigen::VectorXd coeffs,
        Eigen::VectorXd u,
        double variance,
        Eigen::VectorXd marker_variance,
        std::optional<MarkerAllocation> group,
        std::optional<Assignment> sign)
        : type(type),
          coeffs(std::move(coeffs)),
          u(std::move(u)),
          variance(variance),
          marker_variance(std::move(marker_variance)),
          group(std::move(group)),
          sign(std::move(sign))
    {
    }

    GeneticState(
        const GeneticEffect& effect,
        const OldGeneticPrior& prior,
        GeneticMode mode);

    GeneticMode type;
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

namespace gelex::bayes::vi
{

struct GeneticState
{
    GeneticState(
        const bayes::GeneticEffect& effect,
        const bayes::OldGeneticPrior& prior,
        GeneticMode mode);

    GeneticMode type;
    Eigen::VectorXd coeffs;
    Eigen::VectorXd sigma2;
    Eigen::VectorXd u;
    double variance{};
    double heritability{};
    Eigen::VectorXd marker_variance;
};

}  // namespace gelex::bayes::vi

#endif  // GELEX_ALGO_INFER_MCMC_STATE_H_
