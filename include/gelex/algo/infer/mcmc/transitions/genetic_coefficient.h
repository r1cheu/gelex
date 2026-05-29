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

#ifndef GELEX_ALGO_INFER_MCMC_TRANSITIONS_GENETIC_COEFFICIENT_H_
#define GELEX_ALGO_INFER_MCMC_TRANSITIONS_GENETIC_COEFFICIENT_H_

#include <random>

#include <Eigen/Core>

#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class SingleSharedGaussianTransition
{
   public:
    SingleSharedGaussianTransition(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual);

    auto prepare() -> void;

    auto update(
        Eigen::Index marker_index,
        Eigen::Ref<const Eigen::VectorXd> column,
        double xtx_diag_i,
        std::mt19937_64& rng) -> void;

    auto finish() -> void;

   private:
    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;
    double& variance_;
    stats::NormalSampler<double> normal_;
};

class SinglePerMarkerGaussianTransition
{
   public:
    SinglePerMarkerGaussianTransition(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual);

    auto prepare() -> void;

    auto update(
        Eigen::Index marker_index,
        Eigen::Ref<const Eigen::VectorXd> column,
        double xtx_diag_i,
        std::mt19937_64& rng) -> void;

    auto finish() -> void;

   private:
    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;
    Eigen::VectorXd& variance_;
    stats::NormalSampler<double> normal_;
};

class SingleSharedSpikeSlabTransition
{
   public:
    SingleSharedSpikeSlabTransition(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual);

    auto prepare() -> void;

    auto update(
        Eigen::Index marker_index,
        Eigen::Ref<const Eigen::VectorXd> column,
        double xtx_diag_i,
        std::mt19937_64& rng) -> void;

    auto finish() -> void;

   private:
    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;
    double& variance_;
    bayes::ProportionState& proportion_;
    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
};

class SinglePerMarkerSpikeSlabTransition
{
   public:
    SinglePerMarkerSpikeSlabTransition(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual);

    auto prepare() -> void;

    auto update(
        Eigen::Index marker_index,
        Eigen::Ref<const Eigen::VectorXd> column,
        double xtx_diag_i,
        std::mt19937_64& rng) -> void;

    auto finish() -> void;

   private:
    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;
    Eigen::VectorXd& variance_;
    bayes::ProportionState& proportion_;
    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
};

class SingleScaledMixtureTransition
{
   public:
    SingleScaledMixtureTransition(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual);

    auto prepare() -> void;

    auto update(
        Eigen::Index marker_index,
        Eigen::Ref<const Eigen::VectorXd> column,
        double xtx_diag_i,
        std::mt19937_64& rng) -> void;

    auto finish() -> void;

   private:
    static constexpr int kMaxMixtureComponents = 5;

    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;
    double& variance_;
    bayes::ProportionState& proportion_;
    bayes::ComponentState& component_;
    Eigen::VectorXd multiplier_;
    Eigen::VectorXd marker_variances_;
    Eigen::VectorXd logpi_;
    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::Array<double, kMaxMixtureComponents, 1> scale_means_;
    Eigen::Array<double, kMaxMixtureComponents, 1> scale_vars_;
    Eigen::Array<double, kMaxMixtureComponents, 1> scale_log_likelihoods_;
};

class JointGaussianMixtureTransition
{
   public:
    JointGaussianMixtureTransition(
        const bayes::JointGeneticPrior& prior,
        bayes::JointGeneticBlockState& block,
        bayes::ResidualState& residual);

    auto prepare() -> void;

    auto update(
        Eigen::Index marker_index,
        Eigen::Ref<const Eigen::VectorXd> additive_column,
        double additive_xtx_diag_i,
        Eigen::Ref<const Eigen::VectorXd> dominance_column,
        double dominance_xtx_diag_i,
        std::mt19937_64& rng) -> void;

    auto finish() -> void;

   private:
    bayes::GeneticState& additive_;
    bayes::GeneticState& dominance_;
    bayes::ResidualState& residual_;
    bayes::JointSharedVarianceStateCap& variance_;
    bayes::ProportionState& proportion_;
    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_TRANSITIONS_GENETIC_COEFFICIENT_H_
