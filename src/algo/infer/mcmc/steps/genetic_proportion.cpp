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

#include "gelex/algo/infer/mcmc/steps/genetic_proportion.h"
#include <random>
#include <type_traits>
#include <variant>

#include "gelex/exception.h"
#include "gelex/model/bayes/state.h"

namespace gelex::mcmc
{

SingleProportionStep::SingleProportionStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : proportion_(
          std::visit(
              [](auto& state) -> bayes::SampledProportionState&
              {
                  using State = std::decay_t<decltype(state)>;
                  if constexpr (
                      std::is_same_v<
                          State,
                          bayes::SingleSampledSharedSpikeSlabGaussianState>
                      || std::is_same_v<
                          State,
                          bayes::SingleSampledPerMarkerSpikeSlabGaussianState>
                      || std::is_same_v<
                          State,
                          bayes::SingleSampledScaledMixtureGaussianState>)
                  {
                      return state.proportion();
                  }
                  else
                  {
                      throw GelexException(
                          "SingleProportionStep requires sampled proportion "
                          "state");
                  }
              },
              block.prior_state())),
      dirichlet_(
          std::visit(
              [](const auto& value) -> const Eigen::VectorXd&
              {
                  using Prior = std::decay_t<decltype(value)>;
                  if constexpr (
                      std::is_same_v<
                          Prior,
                          bayes::SingleSampledSharedSpikeSlabGaussianPrior>
                      || std::is_same_v<
                          Prior,
                          bayes::SingleSampledPerMarkerSpikeSlabGaussianPrior>
                      || std::is_same_v<
                          Prior,
                          bayes::SingleSampledScaledMixtureGaussianPrior>)
                  {
                      return value.proportion()
                          .parameter()
                          .prior()
                          .concentration();
                  }
                  else
                  {
                      throw GelexException(
                          "SingleProportionStep requires sampled proportion "
                          "prior");
                  }
              },
              prior)),
      rng_(rng)
{
}

auto SingleProportionStep::step() -> void
{
    dirichlet_.reset();
    proportion_.value = dirichlet_(proportion_.assignment.count, rng_);
}

JointProportionStep::JointProportionStep(
    bayes::JointGeneticBlockState& block,
    const bayes::JointGeneticPrior& prior,
    std::mt19937_64& rng)
    : proportion_(
          std::visit(
              [](auto& state) -> bayes::SampledProportionState&
              {
                  using State = std::decay_t<decltype(state)>;
                  if constexpr (
                      std::is_same_v<
                          State,
                          bayes::JointSampledGaussianMixtureState>)
                  {
                      return state.proportion();
                  }
                  else
                  {
                      throw GelexException(
                          "JointProportionStep requires sampled proportion "
                          "state");
                  }
              },
              block.prior_state())),
      dirichlet_(
          std::visit(
              [](const auto& value) -> const Eigen::VectorXd&
              {
                  using Prior = std::decay_t<decltype(value)>;
                  if constexpr (
                      std::is_same_v<
                          Prior,
                          bayes::JointSampledGaussianMixturePrior>)
                  {
                      return value.proportion()
                          .parameter()
                          .prior()
                          .concentration();
                  }
                  else
                  {
                      throw GelexException(
                          "JointProportionStep requires sampled proportion "
                          "prior");
                  }
              },
              prior)),
      rng_(rng)
{
}

auto JointProportionStep::step() -> void
{
    dirichlet_.reset();
    proportion_.value = dirichlet_(proportion_.assignment.count, rng_);
}

}  // namespace gelex::mcmc
