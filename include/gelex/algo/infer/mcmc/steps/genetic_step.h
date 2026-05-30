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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_STEP_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_STEP_H_

#include <array>
#include <optional>
#include <random>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/prior_state_values.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/stats/conjugate_prior.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class SingleSharedGaussianStep final : public Step
{
   public:
    SingleSharedGaussianStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    const bayes::GeneticDesign& design_;

    stats::ScaledInvChi2Sampler<double> variance_sampler_;
    double& variance_;

    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;

    stats::NormalSampler<double> normal_;

    std::mt19937_64& rng_;
};

class SinglePerMarkerGaussianStep final : public Step
{
   public:
    SinglePerMarkerGaussianStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    const bayes::GeneticDesign& design_;

    stats::ScaledInvChi2Sampler<double> variance_sampler_;
    Eigen::VectorXd& variance_;

    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;

    stats::NormalSampler<double> normal_;

    std::mt19937_64& rng_;
};

class SingleSharedSpikeSlabStep final : public Step
{
   public:
    SingleSharedSpikeSlabStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    const bayes::GeneticDesign& design_;

    stats::ScaledInvChi2Sampler<double> variance_sampler_;
    double& variance_;
    Eigen::VectorXi& assignment_;
    Eigen::VectorXd& proportion_;
    std::optional<stats::DirichletSampler<double>> proportion_sampler_;

    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;

    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
    Eigen::VectorXi proportion_count_;

    std::mt19937_64& rng_;
};

class SinglePerMarkerSpikeSlabStep final : public Step
{
   public:
    SinglePerMarkerSpikeSlabStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    const bayes::GeneticDesign& design_;

    stats::ScaledInvChi2Sampler<double> variance_sampler_;
    Eigen::VectorXd& variance_;
    Eigen::VectorXi& assignment_;
    Eigen::VectorXd& proportion_;
    std::optional<stats::DirichletSampler<double>> proportion_sampler_;

    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;

    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
    Eigen::VectorXi proportion_count_;

    std::mt19937_64& rng_;
};

class SingleScaledMixtureStep final : public Step
{
   public:
    SingleScaledMixtureStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    static constexpr int kMaxMixtureComponents = 5;

    const bayes::GeneticDesign& design_;

    stats::ScaledInvChi2Sampler<double> variance_sampler_;
    double& variance_;
    Eigen::VectorXi& assignment_;
    Eigen::VectorXd& proportion_;
    std::optional<stats::DirichletSampler<double>> proportion_sampler_;
    bayes::ComponentState& component_;
    Eigen::VectorXd multiplier_;

    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;

    Eigen::VectorXd marker_variances_;
    Eigen::VectorXd logpi_;
    Eigen::VectorXi proportion_count_;
    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::Array<double, kMaxMixtureComponents, 1> scale_means_;
    Eigen::Array<double, kMaxMixtureComponents, 1> scale_vars_;
    Eigen::Array<double, kMaxMixtureComponents, 1> scale_log_likelihoods_;

    std::mt19937_64& rng_;
};

class JointGaussianMixtureStep final : public Step
{
   public:
    JointGaussianMixtureStep(
        const bayes::GeneticDesign& additive,
        const bayes::GeneticDesign& dominance,
        const bayes::JointGeneticPrior& prior,
        bayes::JointGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    const bayes::GeneticDesign& additive_design_;
    const bayes::GeneticDesign& dominance_design_;

    std::array<stats::ScaledInvChi2Sampler<double>, 2> variance_samplers_;
    std::array<double*, 2> variance_;
    Eigen::VectorXi& assignment_;
    Eigen::VectorXd& proportion_;
    std::optional<stats::DirichletSampler<double>> proportion_sampler_;

    bayes::GeneticState& additive_;
    bayes::GeneticState& dominance_;
    bayes::ResidualState& residual_;

    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
    Eigen::VectorXi proportion_count_;

    std::mt19937_64& rng_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_STEP_H_
