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
#include <variant>

#include <Eigen/Core>

#include "gelex/bayes/state.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

SingleSharedVarStep::SingleSharedVarStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          std::get<bayes::SingleSharedGaussianPrior>(prior)
              .variance()
              .parameter()
              .prior()),
      variance_(
          std::get<bayes::SingleSharedGaussianState>(block.prior_state())
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

SingleSharedSpikeSlabVarStep::SingleSharedSpikeSlabVarStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          std::get<bayes::SingleSharedSpikeSlabGaussianPrior>(prior)
              .variance()
              .parameter()
              .prior()),
      variance_(
          std::get<bayes::SingleSharedSpikeSlabGaussianState>(
              block.prior_state())
              .variance()),
      assignment_(
          std::get<bayes::SingleSharedSpikeSlabGaussianState>(
              block.prior_state())
              .assignment()),
      rng_(rng)
{
}

auto SingleSharedSpikeSlabVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    Eigen::Index n = 0;
    double sum_squares = 0.0;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const int component = assignment_(i);
        if (component == 0)
        {
            continue;
        }
        ++n;
        const double coeff = coeffs(i);
        sum_squares += coeff * coeff;
    }
    variance_ = sampler_({n, sum_squares}, rng_);
}

SingleScaledMixtureVarStep::SingleScaledMixtureVarStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          std::get<bayes::SingleScaledMixtureGaussianPrior>(prior)
              .variance()
              .parameter()
              .prior()),
      variance_(
          std::get<bayes::SingleScaledMixtureGaussianState>(block.prior_state())
              .variance()),
      assignment_(
          std::get<bayes::SingleScaledMixtureGaussianState>(block.prior_state())
              .assignment()),
      multiplier_(
          std::get<bayes::SingleScaledMixtureGaussianPrior>(prior)
              .multiplier()),
      rng_(rng)
{
}

auto SingleScaledMixtureVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    Eigen::Index n = 0;
    double sum_squares = 0.0;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const int component = assignment_(i);
        if (component == 0)
        {
            continue;
        }
        ++n;
        const double coeff = coeffs(i);
        sum_squares += (coeff * coeff) / multiplier_(component);
    }
    variance_ = sampler_({n, sum_squares}, rng_);
}

SinglePerMarkerVarStep::SinglePerMarkerVarStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          std::get<bayes::SinglePerMarkerGaussianPrior>(prior)
              .variance()
              .parameter()
              .prior()),
      variance_(
          std::get<bayes::SinglePerMarkerGaussianState>(block.prior_state())
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

SinglePerMarkerSpikeSlabVarStep::SinglePerMarkerSpikeSlabVarStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          std::get<bayes::SinglePerMarkerSpikeSlabGaussianPrior>(prior)
              .variance()
              .parameter()
              .prior()),
      variance_(
          std::get<bayes::SinglePerMarkerSpikeSlabGaussianState>(
              block.prior_state())
              .variance()),
      assignment_(
          std::get<bayes::SinglePerMarkerSpikeSlabGaussianState>(
              block.prior_state())
              .assignment()),
      rng_(rng)
{
}

auto SinglePerMarkerSpikeSlabVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        if (assignment_(i) == 0)
        {
            continue;
        }
        const double coeff = coeffs(i);
        const double sum_squares = coeff * coeff;
        variance_(i) = sampler_({1, sum_squares}, rng_);
    }
}

JointMixtureVarStep::JointMixtureVarStep(
    bayes::JointGeneticBlockState& block,
    const bayes::JointGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      samplers_(
          [&prior]
          {
              const auto& concrete
                  = std::get<bayes::JointGaussianMixturePrior>(prior);
              return std::array{
                  stats::ScaledInvChi2Sampler<double>{
                      concrete.variance(GeneticMode::A).parameter().prior()},
                  stats::ScaledInvChi2Sampler<double>{
                      concrete.variance(GeneticMode::D).parameter().prior()}};
          }()),
      variance_{
          &std::get<bayes::JointGaussianMixtureState>(block.prior_state())
               .variance(GeneticMode::A),
          &std::get<bayes::JointGaussianMixtureState>(block.prior_state())
               .variance(GeneticMode::D)},
      assignment_(
          std::get<bayes::JointGaussianMixtureState>(block.prior_state())
              .assignment()),
      rng_(rng)
{
}

auto JointMixtureVarStep::step() -> void
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (auto mode : modes)
    {
        auto& sampler = samplers_[std::to_underlying(mode)];
        sampler.reset();
        const auto& coeffs = block_.state(mode).coeffs;
        Eigen::Index n = 0;
        double sum_squares = 0.0;
        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
        {
            const int component = assignment_(i);
            const bool active = mode == GeneticMode::A
                                    ? (component == 1 || component == 3)
                                    : (component == 2 || component == 3);
            if (!active)
            {
                continue;
            }
            ++n;
            const double coeff = coeffs(i);
            sum_squares += coeff * coeff;
        }
        *variance_[std::to_underlying(mode)] = sampler({n, sum_squares}, rng_);
    }
}

}  // namespace gelex::mcmc
