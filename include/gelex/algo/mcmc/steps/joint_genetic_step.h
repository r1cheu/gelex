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

#ifndef GELEX_ALGO_MCMC_STEPS_JOINT_GENETIC_STEP_H_
#define GELEX_ALGO_MCMC_STEPS_JOINT_GENETIC_STEP_H_

#include <Eigen/Core>
#include <array>
#include <cstdint>
#include <optional>
#include <random>
#include <vector>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/half_normal_prior.h"
#include "gelex/bayes/genetic/half_normal_prior_state.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/stats/beta_sampler.h"
#include "gelex/infra/stats/dirichlet_sampler.h"
#include "gelex/infra/stats/half_normal_sampler.h"
#include "gelex/infra/stats/normal_sampler.h"
#include "gelex/infra/stats/scaled_inv_chi2_sampler.h"
#include "gelex/types/categorical_vector.h"

namespace gelex
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class JointGaussianMixtureStep final
{
   public:
    JointGaussianMixtureStep(
        const bayes::GeneticDesign& additive,
        const bayes::GeneticDesign& dominance,
        const bayes::JointGeneticPrior& prior,
        bayes::JointGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void;

   private:
    JointGaussianMixtureStep(
        const bayes::GeneticDesign& additive,
        const bayes::GeneticDesign& dominance,
        const bayes::JointGaussianMixturePrior& prior,
        bayes::JointGaussianMixtureState& prior_state,
        bayes::JointGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    const bayes::GeneticDesign& additive_design_;
    const bayes::GeneticDesign& dominance_design_;
    std::vector<int64_t> valid_indices_;

    std::array<ScaledInvChi2Sampler<double>, 2> variance_samplers_;
    std::array<double*, 2> variance_;
    CategoricalVector& assignment_;
    Eigen::VectorXd& proportion_;
    std::optional<DirichletSampler<double>> proportion_sampler_;

    bayes::GeneticState& additive_;
    bayes::GeneticState& dominance_;
    bayes::ResidualState& residual_;

    NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
    Eigen::VectorXi proportion_count_;

    std::mt19937_64& rng_;
};

class JointHalfNormalMixtureStep final
{
   public:
    JointHalfNormalMixtureStep(
        const bayes::GeneticDesign& additive,
        const bayes::GeneticDesign& dominance,
        const bayes::JointGeneticPrior& prior,
        bayes::JointGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    auto step() -> void;

   private:
    JointHalfNormalMixtureStep(
        const bayes::GeneticDesign& additive,
        const bayes::GeneticDesign& dominance,
        const bayes::JointHalfNormalMixturePrior& prior,
        bayes::JointHalfNormalMixtureState& prior_state,
        bayes::JointGeneticBlockState& block,
        bayes::ResidualState& residual,
        std::mt19937_64& rng);

    const bayes::GeneticDesign& additive_design_;
    const bayes::GeneticDesign& dominance_design_;
    std::vector<int64_t> valid_indices_;

    std::array<ScaledInvChi2Sampler<double>, 2> variance_samplers_;
    std::array<double*, 2> variance_;
    CategoricalVector& assignment_;
    Eigen::VectorXd& proportion_;
    std::optional<DirichletSampler<double>> proportion_sampler_;
    bayes::DominanceSignState& dominance_sign_;

    bayes::GeneticState& additive_;
    bayes::GeneticState& dominance_;
    bayes::ResidualState& residual_;

    NormalSampler<double> normal_;
    HalfNormalSampler<double> half_normal_;
    BetaSampler<double> sign_sampler_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
    Eigen::VectorXi proportion_count_;

    std::mt19937_64& rng_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex

#endif  // GELEX_ALGO_MCMC_STEPS_JOINT_GENETIC_STEP_H_
