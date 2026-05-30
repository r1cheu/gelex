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

SingleFixedSharedSpikeSlabVarStep::SingleFixedSharedSpikeSlabVarStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          std::get<bayes::SingleFixedSharedSpikeSlabGaussianPrior>(prior)
              .variance()
              .parameter()
              .prior()),
      variance_(
          std::get<bayes::SingleFixedSharedSpikeSlabGaussianState>(
              block.prior_state())
              .variance()),
      assignment_(
          std::get<bayes::SingleFixedSharedSpikeSlabGaussianState>(
              block.prior_state())
              .assignment()),
      rng_(rng)
{
}

auto SingleFixedSharedSpikeSlabVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    const Eigen::Index n
        = assignment_.count.tail(assignment_.count.size() - 1).sum();
    double sum_squares = 0.0;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const int component = assignment_.assignment(i);
        if (component == 0)
        {
            continue;
        }
        const double coeff = coeffs(i);
        sum_squares += coeff * coeff;
    }
    variance_ = sampler_({n, sum_squares}, rng_);
}

SingleSampledSharedSpikeSlabVarStep::SingleSampledSharedSpikeSlabVarStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          std::get<bayes::SingleSampledSharedSpikeSlabGaussianPrior>(prior)
              .variance()
              .parameter()
              .prior()),
      variance_(
          std::get<bayes::SingleSampledSharedSpikeSlabGaussianState>(
              block.prior_state())
              .variance()),
      assignment_(
          std::get<bayes::SingleSampledSharedSpikeSlabGaussianState>(
              block.prior_state())
              .assignment()),
      rng_(rng)
{
}

auto SingleSampledSharedSpikeSlabVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    const Eigen::Index n
        = assignment_.count.tail(assignment_.count.size() - 1).sum();
    double sum_squares = 0.0;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const int component = assignment_.assignment(i);
        if (component == 0)
        {
            continue;
        }
        const double coeff = coeffs(i);
        sum_squares += coeff * coeff;
    }
    variance_ = sampler_({n, sum_squares}, rng_);
}

SingleFixedScaledMixtureVarStep::SingleFixedScaledMixtureVarStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          std::get<bayes::SingleFixedScaledMixtureGaussianPrior>(prior)
              .variance()
              .parameter()
              .prior()),
      variance_(
          std::get<bayes::SingleFixedScaledMixtureGaussianState>(
              block.prior_state())
              .variance()),
      assignment_(
          std::get<bayes::SingleFixedScaledMixtureGaussianState>(
              block.prior_state())
              .assignment()),
      multiplier_(
          std::get<bayes::SingleFixedScaledMixtureGaussianPrior>(prior)
              .multiplier()),
      rng_(rng)
{
}

auto SingleFixedScaledMixtureVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    const Eigen::Index n
        = assignment_.count.tail(assignment_.count.size() - 1).sum();
    double sum_squares = 0.0;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const int component = assignment_.assignment(i);
        if (component == 0)
        {
            continue;
        }
        const double coeff = coeffs(i);
        sum_squares += (coeff * coeff) / multiplier_(component);
    }
    variance_ = sampler_({n, sum_squares}, rng_);
}

SingleSampledScaledMixtureVarStep::SingleSampledScaledMixtureVarStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          std::get<bayes::SingleSampledScaledMixtureGaussianPrior>(prior)
              .variance()
              .parameter()
              .prior()),
      variance_(
          std::get<bayes::SingleSampledScaledMixtureGaussianState>(
              block.prior_state())
              .variance()),
      assignment_(
          std::get<bayes::SingleSampledScaledMixtureGaussianState>(
              block.prior_state())
              .assignment()),
      multiplier_(
          std::get<bayes::SingleSampledScaledMixtureGaussianPrior>(prior)
              .multiplier()),
      rng_(rng)
{
}

auto SingleSampledScaledMixtureVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    const Eigen::Index n
        = assignment_.count.tail(assignment_.count.size() - 1).sum();
    double sum_squares = 0.0;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const int component = assignment_.assignment(i);
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

SingleFixedPerMarkerSpikeSlabVarStep::SingleFixedPerMarkerSpikeSlabVarStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          std::get<bayes::SingleFixedPerMarkerSpikeSlabGaussianPrior>(prior)
              .variance()
              .parameter()
              .prior()),
      variance_(
          std::get<bayes::SingleFixedPerMarkerSpikeSlabGaussianState>(
              block.prior_state())
              .variance()),
      assignment_(
          std::get<bayes::SingleFixedPerMarkerSpikeSlabGaussianState>(
              block.prior_state())
              .assignment()),
      rng_(rng)
{
}

auto SingleFixedPerMarkerSpikeSlabVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        if (assignment_.assignment(i) == 0)
        {
            continue;
        }
        const double coeff = coeffs(i);
        const double sum_squares = coeff * coeff;
        variance_(i) = sampler_({1, sum_squares}, rng_);
    }
}

SingleSampledPerMarkerSpikeSlabVarStep::SingleSampledPerMarkerSpikeSlabVarStep(
    bayes::SingleGeneticBlockState& block,
    const bayes::SingleGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      sampler_(
          std::get<bayes::SingleSampledPerMarkerSpikeSlabGaussianPrior>(prior)
              .variance()
              .parameter()
              .prior()),
      variance_(
          std::get<bayes::SingleSampledPerMarkerSpikeSlabGaussianState>(
              block.prior_state())
              .variance()),
      assignment_(
          std::get<bayes::SingleSampledPerMarkerSpikeSlabGaussianState>(
              block.prior_state())
              .assignment()),
      rng_(rng)
{
}

auto SingleSampledPerMarkerSpikeSlabVarStep::step() -> void
{
    sampler_.reset();
    const auto& coeffs = block_.state().coeffs;
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        if (assignment_.assignment(i) == 0)
        {
            continue;
        }
        const double coeff = coeffs(i);
        const double sum_squares = coeff * coeff;
        variance_(i) = sampler_({1, sum_squares}, rng_);
    }
}

JointFixedMixtureVarStep::JointFixedMixtureVarStep(
    bayes::JointGeneticBlockState& block,
    const bayes::JointGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      samplers_(
          [&prior]
          {
              const auto& concrete
                  = std::get<bayes::JointFixedGaussianMixturePrior>(prior);
              return std::array{
                  stats::ScaledInvChi2Sampler<double>{
                      concrete.variance(GeneticMode::A).parameter().prior()},
                  stats::ScaledInvChi2Sampler<double>{
                      concrete.variance(GeneticMode::D).parameter().prior()}};
          }()),
      variance_{
          &std::get<bayes::JointFixedGaussianMixtureState>(block.prior_state())
               .variance(GeneticMode::A),
          &std::get<bayes::JointFixedGaussianMixtureState>(block.prior_state())
               .variance(GeneticMode::D)},
      assignment_(
          std::get<bayes::JointFixedGaussianMixtureState>(block.prior_state())
              .assignment()),
      rng_(rng)
{
}

auto JointFixedMixtureVarStep::step() -> void
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (auto mode : modes)
    {
        auto& sampler = samplers_[std::to_underlying(mode)];
        sampler.reset();
        const auto& coeffs = block_.state(mode).coeffs;
        const Eigen::Index n
            = mode == GeneticMode::A
                  ? assignment_.count(1) + assignment_.count(3)
                  : assignment_.count(2) + assignment_.count(3);
        double sum_squares = 0.0;
        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
        {
            const int component = assignment_.assignment(i);
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
        *variance_[std::to_underlying(mode)] = sampler({n, sum_squares}, rng_);
    }
}

JointSampledMixtureVarStep::JointSampledMixtureVarStep(
    bayes::JointGeneticBlockState& block,
    const bayes::JointGeneticPrior& prior,
    std::mt19937_64& rng)
    : block_(block),
      samplers_(
          [&prior]
          {
              const auto& concrete
                  = std::get<bayes::JointSampledGaussianMixturePrior>(prior);
              return std::array{
                  stats::ScaledInvChi2Sampler<double>{
                      concrete.variance(GeneticMode::A).parameter().prior()},
                  stats::ScaledInvChi2Sampler<double>{
                      concrete.variance(GeneticMode::D).parameter().prior()}};
          }()),
      variance_{
          &std::get<bayes::JointSampledGaussianMixtureState>(
               block.prior_state())
               .variance(GeneticMode::A),
          &std::get<bayes::JointSampledGaussianMixtureState>(
               block.prior_state())
               .variance(GeneticMode::D)},
      assignment_(
          std::get<bayes::JointSampledGaussianMixtureState>(block.prior_state())
              .assignment()),
      rng_(rng)
{
}

auto JointSampledMixtureVarStep::step() -> void
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (auto mode : modes)
    {
        auto& sampler = samplers_[std::to_underlying(mode)];
        sampler.reset();
        const auto& coeffs = block_.state(mode).coeffs;
        const Eigen::Index n
            = mode == GeneticMode::A
                  ? assignment_.count(1) + assignment_.count(3)
                  : assignment_.count(2) + assignment_.count(3);
        double sum_squares = 0.0;
        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
        {
            const int component = assignment_.assignment(i);
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
        *variance_[std::to_underlying(mode)] = sampler({n, sum_squares}, rng_);
    }
}

}  // namespace gelex::mcmc
