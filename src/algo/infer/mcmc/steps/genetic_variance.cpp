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

#include "gelex/algo/infer/mcmc/steps/genetic_variance.h"

#include <array>
#include <random>
#include <utility>

#include <Eigen/Core>

#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

SingleSharedVarStep::SingleSharedVarStep(
    const bayes::SingleGeneticPrior& prior,
    bayes::SingleGeneticBlockState& block,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          prior.get<bayes::SingleSharedMarkerVarianceCap>()
              .variance()
              .parameter()
              .prior()
              .degrees_of_freedom(),
          prior.get<bayes::SingleSharedMarkerVarianceCap>()
              .variance()
              .parameter()
              .prior()
              .scale()),
      variance_(block.prior_state()
                    .get<bayes::SingleSharedVarianceStateCap>()
                    .variance()),
      rng_(rng)
{
}

auto SingleSharedVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    variance_ = sampler_({coeffs.size(), coeffs.squaredNorm()}, rng_);
}

SingleSharedMixtureVarStep::SingleSharedMixtureVarStep(
    const bayes::SingleGeneticPrior& prior,
    bayes::SingleGeneticBlockState& block,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          prior.get<bayes::SingleSharedMarkerVarianceCap>()
              .variance()
              .parameter()
              .prior()
              .degrees_of_freedom(),
          prior.get<bayes::SingleSharedMarkerVarianceCap>()
              .variance()
              .parameter()
              .prior()
              .scale()),
      variance_(block.prior_state()
                    .get<bayes::SingleSharedVarianceStateCap>()
                    .variance()),
      proportion_(block.prior_state()
                      .get<bayes::SingleProportionStateCap>()
                      .proportion()),
      rng_(rng)
{
}

auto SingleSharedMixtureVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    const Eigen::Index n
        = proportion_.count.tail(proportion_.count.size() - 1).sum();
    double sum_squares = 0.0;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const int component = proportion_.assignment(i);
        if (component == 0)
        {
            continue;
        }
        const double coeff = coeffs(i);
        sum_squares += coeff * coeff;
    }
    variance_ = sampler_({n, sum_squares}, rng_);
}

SingleSharedScaledMixtureVarStep::SingleSharedScaledMixtureVarStep(
    const bayes::SingleGeneticPrior& prior,
    bayes::SingleGeneticBlockState& block,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          prior.get<bayes::SingleSharedMarkerVarianceCap>()
              .variance()
              .parameter()
              .prior()
              .degrees_of_freedom(),
          prior.get<bayes::SingleSharedMarkerVarianceCap>()
              .variance()
              .parameter()
              .prior()
              .scale()),
      variance_(block.prior_state()
                    .get<bayes::SingleSharedVarianceStateCap>()
                    .variance()),
      proportion_(block.prior_state()
                      .get<bayes::SingleProportionStateCap>()
                      .proportion()),
      multiplier_(prior.get<bayes::SingleMultiplierCap>().multiplier()),
      rng_(rng)
{
}

auto SingleSharedScaledMixtureVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    const Eigen::Index n
        = proportion_.count.tail(proportion_.count.size() - 1).sum();
    double sum_squares = 0.0;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const int component = proportion_.assignment(i);
        if (component == 0)
        {
            continue;
        }
        const double coeff = coeffs(i);
        sum_squares += (coeff * coeff) / multiplier_(component);
    }
    variance_ = sampler_({n, sum_squares}, rng_);
}

SinglePerMarkerVarStep::SinglePerMarkerVarStep(
    const bayes::SingleGeneticPrior& prior,
    bayes::SingleGeneticBlockState& block,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          prior.get<bayes::SinglePerMarkerVarianceCap>()
              .variance()
              .parameter()
              .prior()
              .degrees_of_freedom(),
          prior.get<bayes::SinglePerMarkerVarianceCap>()
              .variance()
              .parameter()
              .prior()
              .scale()),
      variance_(block.prior_state()
                    .get<bayes::SinglePerMarkerVarianceStateCap>()
                    .variance()),
      rng_(rng)
{
}

auto SinglePerMarkerVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const double coeff = coeffs(i);
        const double sum_squares = coeff * coeff;
        variance_(i) = sampler_({1, sum_squares}, rng_);
    }
}

SinglePerMarkerMixtureVarStep::SinglePerMarkerMixtureVarStep(
    const bayes::SingleGeneticPrior& prior,
    bayes::SingleGeneticBlockState& block,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          prior.get<bayes::SinglePerMarkerVarianceCap>()
              .variance()
              .parameter()
              .prior()
              .degrees_of_freedom(),
          prior.get<bayes::SinglePerMarkerVarianceCap>()
              .variance()
              .parameter()
              .prior()
              .scale()),
      variance_(block.prior_state()
                    .get<bayes::SinglePerMarkerVarianceStateCap>()
                    .variance()),
      proportion_(block.prior_state()
                      .get<bayes::SingleProportionStateCap>()
                      .proportion()),
      rng_(rng)
{
}

auto SinglePerMarkerMixtureVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        if (proportion_.assignment(i) == 0)
        {
            continue;
        }
        const double coeff = coeffs(i);
        const double sum_squares = coeff * coeff;
        variance_(i) = sampler_({1, sum_squares}, rng_);
    }
}

JointSharedMixtureVarStep::JointSharedMixtureVarStep(
    const bayes::JointGeneticPrior& prior,
    bayes::JointGeneticBlockState& block,
    std::mt19937_64& rng)
    : block_(block),
      samplers_(
          [&prior]
          {
              const auto& variance
                  = prior.get<bayes::JointSharedMarkerVarianceCap>();
              const auto& additive
                  = variance.variance(GeneticMode::A).parameter().prior();
              const auto& dominance
                  = variance.variance(GeneticMode::D).parameter().prior();
              return std::array{
                  stats::ScaledInvChi2Sampler<double>{
                      additive.degrees_of_freedom(), additive.scale()},
                  stats::ScaledInvChi2Sampler<double>{
                      dominance.degrees_of_freedom(), dominance.scale()}};
          }()),
      variance_(block.prior_state().get<bayes::JointSharedVarianceStateCap>()),
      proportion_(block.prior_state()
                      .get<bayes::JointProportionStateCap>()
                      .proportion()),
      rng_(rng)
{
}

auto JointSharedMixtureVarStep::step() -> void
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (auto mode : modes)
    {
        auto& sampler = samplers_[std::to_underlying(mode)];
        sampler.reset();
        const auto& coeffs = block_.state(mode).coeffs;
        const Eigen::Index n
            = mode == GeneticMode::A
                  ? proportion_.count(1) + proportion_.count(3)
                  : proportion_.count(2) + proportion_.count(3);
        double sum_squares = 0.0;
        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
        {
            const int component = proportion_.assignment(i);
            const bool active = mode == GeneticMode::A
                                    ? (component == 1 || component == 3)
                                    : (component == 2 || component == 3);
            if (!active)
            {
                continue;
            }
            const double coeff = coeffs(i);
            sum_squares += coeff * coeff;
        }
        variance_.variance(mode) = sampler({n, sum_squares}, rng_);
    }
}

}  // namespace gelex::mcmc
