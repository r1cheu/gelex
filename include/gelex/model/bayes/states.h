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
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/effects.h"
#include "gelex/types/fixed_effects.h"

namespace gelex::bayes
{

using TrackerVector = Eigen::VectorX<int8_t>;

struct Pi
{
    Eigen::VectorXd proportion;
    Eigen::VectorXi count;
};

struct MixtureState
{
    explicit MixtureState(const GeneticEffect& effect)
        : tracker(TrackerVector::Zero(bayes::get_cols(effect.X))),
          pi{effect.mixture->init_proportion,
             Eigen::VectorXi::Zero(effect.mixture->init_proportion.size())}
    {
        if (const auto num_components = effect.mixture->init_proportion.size();
            num_components > 2)
        {
            const auto num_samples = bayes::get_rows(effect.X);
            component_u.resize(num_components - 1);
            for (auto& vec : component_u)
            {
                vec = Eigen::VectorXd::Zero(num_samples);
            }
            component_variance = Eigen::VectorXd::Zero(num_components - 1);
        }
    }

    TrackerVector tracker;
    Pi pi;
    std::vector<Eigen::VectorXd> component_u;
    Eigen::VectorXd component_variance;
};

struct SignState
{
    explicit SignState(Eigen::Index n_snps)
        : tracker(TrackerVector::Zero(n_snps))
    {
    }
    SignState() = default;

    TrackerVector tracker;
    double positive_prob{0.5};
};

struct FixedState
{
    explicit FixedState(const FixedEffect& effect)
        : coeffs(Eigen::VectorXd::Zero(effect.X.cols())) {};
    Eigen::VectorXd coeffs;
};

struct RandomState
{
    explicit RandomState(const RandomEffect& effect)
        : coeffs(Eigen::VectorXd::Zero(effect.X.cols())),
          variance{effect.init_variance}
    {
    }

    Eigen::VectorXd coeffs;
    double variance{0.0};
};

struct GeneticState
{
    explicit GeneticState(const GeneticEffect& effect)
        : type(effect.type),
          coeffs(Eigen::VectorXd::Zero(bayes::get_cols(effect.X))),
          u(Eigen::VectorXd::Zero(bayes::get_rows(effect.X))),
          marker_variance(
              Eigen::VectorXd::Constant(
                  effect.marker_variance_size,
                  effect.init_marker_variance))
    {
        if (effect.mixture)
        {
            mixture.emplace(effect);
        }
        if (effect.sign)
        {
            sign.emplace(bayes::get_cols(effect.X));
            sign->positive_prob = effect.sign->init_positive_prob;
        }
    }
    GeneticEffectType type;
    Eigen::VectorXd coeffs;
    Eigen::VectorXd u;

    double variance{};
    double heritability{};
    Eigen::VectorXd marker_variance;

    std::optional<MixtureState> mixture;
    std::optional<SignState> sign;
};

struct ResidualState
{
    Eigen::VectorXd y_adj;
    double variance{0.0};
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_STATES_H_
