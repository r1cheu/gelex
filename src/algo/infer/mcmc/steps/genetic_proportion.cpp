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

#include "gelex/bayes/state.h"
#include "gelex/exception.h"

namespace gelex::mcmc
{

SingleProportionStep::SingleProportionStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : mixture_(
          std::visit(
              [](auto& state) -> bayes::MixtureState&
              {
                  using State = std::decay_t<decltype(state)>;
                  if constexpr (
                      std::is_same_v<
                          State,
                          bayes::SingleSharedSpikeSlabGaussianState>
                      || std::is_same_v<
                          State,
                          bayes::SinglePerMarkerSpikeSlabGaussianState>
                      || std::is_same_v<
                          State,
                          bayes::SingleScaledMixtureGaussianState>)
                  {
                      return state.mixture();
                  }
                  else
                  {
                      throw GelexException(
                          "SingleProportionStep requires mixture "
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
                          bayes::SingleSharedSpikeSlabGaussianPrior>
                      || std::is_same_v<
                          Prior,
                          bayes::SinglePerMarkerSpikeSlabGaussianPrior>
                      || std::is_same_v<
                          Prior,
                          bayes::SingleScaledMixtureGaussianPrior>)
                  {
                      if (value.proportion().prior())
                      {
                          return value.proportion().prior()->concentration();
                      }
                      throw GelexException(
                          "SingleProportionStep requires sampled proportion "
                          "prior");
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
    Eigen::VectorXi count = Eigen::VectorXi::Zero(mixture_.proportion.size());
    for (Eigen::Index i = 0; i < mixture_.assignment.assignment.size(); ++i)
    {
        ++count(mixture_.assignment.assignment(i));
    }
    mixture_.proportion = dirichlet_(count, rng_);
}

JointProportionStep::JointProportionStep(
    bayes::JointGeneticBlockState& block,
    const bayes::JointGeneticPrior& prior,
    std::mt19937_64& rng)
    : mixture_(
          std::visit(
              [](auto& state) -> bayes::MixtureState&
              {
                  using State = std::decay_t<decltype(state)>;
                  if constexpr (
                      std::is_same_v<State, bayes::JointGaussianMixtureState>)
                  {
                      return state.mixture();
                  }
                  else
                  {
                      throw GelexException(
                          "JointProportionStep requires mixture "
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
                      std::is_same_v<Prior, bayes::JointGaussianMixturePrior>)
                  {
                      if (value.proportion().prior())
                      {
                          return value.proportion().prior()->concentration();
                      }
                      throw GelexException(
                          "JointProportionStep requires sampled proportion "
                          "prior");
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
    Eigen::VectorXi count = Eigen::VectorXi::Zero(mixture_.proportion.size());
    for (Eigen::Index i = 0; i < mixture_.assignment.assignment.size(); ++i)
    {
        ++count(mixture_.assignment.assignment(i));
    }
    mixture_.proportion = dirichlet_(count, rng_);
}

}  // namespace gelex::mcmc
