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

#ifndef GELEX_ALGO_MCMC_STEPS_SINGLE_GENETIC_STEP_H_
#define GELEX_ALGO_MCMC_STEPS_SINGLE_GENETIC_STEP_H_

#include <Eigen/Core>
#include <array>
#include <optional>
#include <random>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/legacy_state.h"
#include "gelex/infra/stats/dirichlet_sampler.h"
#include "gelex/infra/stats/normal_sampler.h"
#include "gelex/infra/stats/scaled_inv_chi2_sampler.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class SingleSharedGaussianStep final
{
   public:
    SingleSharedGaussianStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleSharedGaussianPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void;

   private:
    const bayes::GeneticDesign& design_;
    GeneticMode mode_;
    ScaledInvChi2Sampler<double> variance_sampler_;
    bayes::SingleGeneticBlockState& block_;
    bayes::ResidualState& residual_;
    NormalSampler<double> normal_;
    std::mt19937_64& rng_;
};

class SinglePerMarkerGaussianStep final
{
   public:
    SinglePerMarkerGaussianStep(
        const bayes::GeneticDesign& design,
        const bayes::SinglePerMarkerGaussianPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void;

   private:
    const bayes::GeneticDesign& design_;
    GeneticMode mode_;
    ScaledInvChi2Sampler<double> variance_sampler_;
    bayes::SingleGeneticBlockState& block_;
    bayes::ResidualState& residual_;
    NormalSampler<double> normal_;
    std::mt19937_64& rng_;
};

class SingleSharedSpikeSlabStep final
{
   public:
    SingleSharedSpikeSlabStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleSharedSpikeSlabGaussianPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void;

   private:
    const bayes::GeneticDesign& design_;
    GeneticMode mode_;
    ScaledInvChi2Sampler<double> variance_sampler_;
    std::optional<DirichletSampler<double>> proportion_sampler_;
    bayes::SingleGeneticBlockState& block_;
    bayes::ResidualState& residual_;
    NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
    Eigen::VectorXi proportion_count_;
    std::mt19937_64& rng_;
};

class SinglePerMarkerSpikeSlabStep final
{
   public:
    SinglePerMarkerSpikeSlabStep(
        const bayes::GeneticDesign& design,
        const bayes::SinglePerMarkerSpikeSlabGaussianPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void;

   private:
    const bayes::GeneticDesign& design_;
    GeneticMode mode_;
    ScaledInvChi2Sampler<double> variance_sampler_;
    std::optional<DirichletSampler<double>> proportion_sampler_;
    bayes::SingleGeneticBlockState& block_;
    bayes::ResidualState& residual_;
    NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
    Eigen::VectorXi proportion_count_;
    std::mt19937_64& rng_;
};

class SingleScaledMixtureStep final
{
   public:
    SingleScaledMixtureStep(
        const bayes::GeneticDesign& design,
        const bayes::SingleScaledMixtureGaussianPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void;

   private:
    static constexpr int max_mixture_components = 5;

    const bayes::GeneticDesign& design_;
    GeneticMode mode_;
    ScaledInvChi2Sampler<double> variance_sampler_;
    DirichletSampler<double> proportion_sampler_;
    Eigen::VectorXd multiplier_;
    bayes::SingleGeneticBlockState& block_;
    bayes::ResidualState& residual_;
    Eigen::VectorXd marker_variances_;
    Eigen::VectorXd logpi_;
    Eigen::VectorXi proportion_count_;
    NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    std::array<NormalSampler<double>::Posterior, max_mixture_components>
        scale_posts_;
    std::mt19937_64& rng_;
};

// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex

#endif  // GELEX_ALGO_MCMC_STEPS_SINGLE_GENETIC_STEP_H_
