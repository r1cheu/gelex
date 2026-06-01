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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_SINGLE_GENETIC_STEP_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_SINGLE_GENETIC_STEP_H_

#include <array>
#include <optional>
#include <random>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
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
    SingleSharedGaussianStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleSharedGaussianPrior& prior,
        bayes::SingleSharedGaussianState& prior_state,
        bayes::GeneticState& state,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

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
    SinglePerMarkerGaussianStep(
        const bayes::GeneticDesign& design,
        const bayes::SinglePerMarkerGaussianPrior& prior,
        bayes::SinglePerMarkerGaussianState& prior_state,
        bayes::GeneticState& state,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

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
    SingleSharedSpikeSlabStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleSharedSpikeSlabGaussianPrior& prior,
        bayes::SingleSharedSpikeSlabGaussianState& prior_state,
        bayes::GeneticState& state,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

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
    SinglePerMarkerSpikeSlabStep(
        const bayes::GeneticDesign& design,
        const bayes::SinglePerMarkerSpikeSlabGaussianPrior& prior,
        bayes::SinglePerMarkerSpikeSlabGaussianState& prior_state,
        bayes::GeneticState& state,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

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
    static constexpr int MAX_MIXTURE_COMPONENTS = 5;

    SingleScaledMixtureStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleScaledMixtureGaussianPrior& prior,
        bayes::SingleScaledMixtureGaussianState& prior_state,
        bayes::GeneticState& state,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    const bayes::GeneticDesign& design_;

    stats::ScaledInvChi2Sampler<double> variance_sampler_;
    double& variance_;
    Eigen::VectorXi& assignment_;
    Eigen::VectorXd& proportion_;
    stats::DirichletSampler<double> proportion_sampler_;
    bayes::ComponentState& component_;
    Eigen::VectorXd multiplier_;

    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;

    Eigen::VectorXd marker_variances_;
    Eigen::VectorXd logpi_;
    Eigen::VectorXi proportion_count_;
    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    std::array<stats::NormalSampler<double>::Posterior, MAX_MIXTURE_COMPONENTS>
        scale_posts_;
    Eigen::Array<double, MAX_MIXTURE_COMPONENTS, 1> scale_log_likelihoods_;

    std::mt19937_64& rng_;
};

// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_SINGLE_GENETIC_STEP_H_
