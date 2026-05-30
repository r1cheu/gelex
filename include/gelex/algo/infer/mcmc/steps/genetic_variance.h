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
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class SingleSharedVarStep final : public Step
{
   public:
    SingleSharedVarStep(
        bayes::SingleGeneticBlockState& block,
        const bayes::SingleGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    double& variance_;
    std::mt19937_64& rng_;
};

class SingleFixedSharedSpikeSlabVarStep final : public Step
{
   public:
    SingleFixedSharedSpikeSlabVarStep(
        bayes::SingleGeneticBlockState& block,
        const bayes::SingleGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    double& variance_;
    bayes::MixtureAssignmentState& assignment_;
    std::mt19937_64& rng_;
};

class SingleSampledSharedSpikeSlabVarStep final : public Step
{
   public:
    SingleSampledSharedSpikeSlabVarStep(
        bayes::SingleGeneticBlockState& block,
        const bayes::SingleGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    double& variance_;
    bayes::MixtureAssignmentState& assignment_;
    std::mt19937_64& rng_;
};

class SingleFixedScaledMixtureVarStep final : public Step
{
   public:
    SingleFixedScaledMixtureVarStep(
        bayes::SingleGeneticBlockState& block,
        const bayes::SingleGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    double& variance_;
    bayes::MixtureAssignmentState& assignment_;
    const Eigen::VectorXd& multiplier_;
    std::mt19937_64& rng_;
};

class SingleSampledScaledMixtureVarStep final : public Step
{
   public:
    SingleSampledScaledMixtureVarStep(
        bayes::SingleGeneticBlockState& block,
        const bayes::SingleGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    double& variance_;
    bayes::MixtureAssignmentState& assignment_;
    const Eigen::VectorXd& multiplier_;
    std::mt19937_64& rng_;
};

class SinglePerMarkerVarStep final : public Step
{
   public:
    SinglePerMarkerVarStep(
        bayes::SingleGeneticBlockState& block,
        const bayes::SingleGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    Eigen::VectorXd& variance_;
    std::mt19937_64& rng_;
};

class SingleFixedPerMarkerSpikeSlabVarStep final : public Step
{
   public:
    SingleFixedPerMarkerSpikeSlabVarStep(
        bayes::SingleGeneticBlockState& block,
        const bayes::SingleGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    Eigen::VectorXd& variance_;
    bayes::MixtureAssignmentState& assignment_;
    std::mt19937_64& rng_;
};

class SingleSampledPerMarkerSpikeSlabVarStep final : public Step
{
   public:
    SingleSampledPerMarkerSpikeSlabVarStep(
        bayes::SingleGeneticBlockState& block,
        const bayes::SingleGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::SingleGeneticBlockState& block_;
    stats::ScaledInvChi2Sampler<double> sampler_;
    Eigen::VectorXd& variance_;
    bayes::MixtureAssignmentState& assignment_;
    std::mt19937_64& rng_;
};

class JointFixedMixtureVarStep final : public Step
{
   public:
    JointFixedMixtureVarStep(
        bayes::JointGeneticBlockState& block,
        const bayes::JointGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::JointGeneticBlockState& block_;
    std::array<stats::ScaledInvChi2Sampler<double>, 2> samplers_;
    std::array<double*, 2> variance_;
    bayes::MixtureAssignmentState& assignment_;
    std::mt19937_64& rng_;
};

class JointSampledMixtureVarStep final : public Step
{
   public:
    JointSampledMixtureVarStep(
        bayes::JointGeneticBlockState& block,
        const bayes::JointGeneticPrior& prior,
        std::mt19937_64& rng);

    auto step() -> void override;

   private:
    bayes::JointGeneticBlockState& block_;
    std::array<stats::ScaledInvChi2Sampler<double>, 2> samplers_;
    std::array<double*, 2> variance_;
    bayes::MixtureAssignmentState& assignment_;
    std::mt19937_64& rng_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_VARIANCE_H_
