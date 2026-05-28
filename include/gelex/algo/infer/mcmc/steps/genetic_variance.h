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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_VARIANCE_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_VARIANCE_H_

#include <array>
#include <random>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class SingleSharedGeneticVarianceStep final : public Step
{
   public:
    SingleSharedGeneticVarianceStep(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    double& variance_;
    std::mt19937_64& rng_;
};

class SingleSharedMixtureGeneticVarianceStep final : public Step
{
   public:
    SingleSharedMixtureGeneticVarianceStep(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    double& variance_;
    bayes::ProportionState& proportion_;
    std::mt19937_64& rng_;
};

class SingleSharedScaledMixtureGeneticVarianceStep final : public Step
{
   public:
    SingleSharedScaledMixtureGeneticVarianceStep(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    double& variance_;
    bayes::ProportionState& proportion_;
    const Eigen::VectorXd& multiplier_;
    std::mt19937_64& rng_;
};

class SinglePerMarkerGeneticVarianceStep final : public Step
{
   public:
    SinglePerMarkerGeneticVarianceStep(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    Eigen::VectorXd& variance_;
    std::mt19937_64& rng_;
};

class SinglePerMarkerMixtureGeneticVarianceStep final : public Step
{
   public:
    SinglePerMarkerMixtureGeneticVarianceStep(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    Eigen::VectorXd& variance_;
    bayes::ProportionState& proportion_;
    std::mt19937_64& rng_;
};

class JointSharedMixtureGeneticVarianceStep final : public Step
{
   public:
    JointSharedMixtureGeneticVarianceStep(
        const bayes::JointGeneticPrior& prior,
        bayes::JointGeneticBlockState& block,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::JointGeneticBlockState& block_;
    std::array<stats::ScaledInvChi2Sampler<double>, 2> samplers_;
    bayes::JointSharedVarianceStateCap& variance_;
    bayes::ProportionState& proportion_;
    std::mt19937_64& rng_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_VARIANCE_H_
