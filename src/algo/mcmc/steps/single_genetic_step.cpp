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

#include "gelex/algo/mcmc/steps/single_genetic_step.h"

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <optional>
#include <random>
#include <variant>

#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/legacy_genetic_prior.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/infra/stats/dirichlet_sampler.h"
#include "gelex/infra/stats/normal_sampler.h"

#include "algo/mcmc/invariant.h"

namespace gelex
{

SingleSharedGaussianStep::SingleSharedGaussianStep(
    const bayes::GeneticDesign& design,
    const bayes::SingleSharedGaussianPrior& prior,
    bayes::SingleGeneticBlockState& block,
    bayes::ResidualState& residual,
    std::mt19937_64& rng)
    : design_(design),
      mode_(block.mode()),
      variance_sampler_(prior.variance().parameter().prior()),
      block_(block),
      residual_(residual),
      normal_(0.0),
      rng_(rng)
{
    (void)std::get<bayes::SingleSharedGaussianState>(block_.prior_state());
}

auto SingleSharedGaussianStep::step() -> void
{
    auto& prior_state
        = std::get<bayes::SingleSharedGaussianState>(block_.prior_state());
    auto& variance = prior_state.variance();
    auto& state = block_.state();
    const auto& xtx_diag = design_.xtx_diag(mode_);
    auto& coeffs = state.coeffs;

    normal_.reset();
    normal_.set_prior_var(variance);
    variance_sampler_.reset();

    Eigen::Index variance_n = 0;
    double sum_squares = 0.0;
    {
        GeneticSweepAdjustmentGuard sweep{residual_, state};
        for (const Eigen::Index i : design_.valid_indices(mode_))
        {
            const double old_i = coeffs(i);
            GeneticResidualAdjustmentGuard guard{
                design_, mode_, i, state, residual_};
            const double rhs = design_.dot(mode_, i, residual_.y_adj)
                               + (xtx_diag(i) * old_i);
            coeffs(i) = normal_(
                NormalSampler<double>::Kernel{
                    .quadratic = xtx_diag(i),
                    .linear = rhs,
                    .scale = residual_.variance,
                },
                rng_);
            ++variance_n;
            const double coeff = coeffs(i);
            sum_squares += coeff * coeff;
        }
    }
    variance = variance_sampler_({variance_n, sum_squares}, rng_);
    state.variance = detail::vecvar(state.u, detail::VarNormType::Population);
}

SinglePerMarkerGaussianStep::SinglePerMarkerGaussianStep(
    const bayes::GeneticDesign& design,
    const bayes::SinglePerMarkerGaussianPrior& prior,
    bayes::SingleGeneticBlockState& block,
    bayes::ResidualState& residual,
    std::mt19937_64& rng)
    : design_(design),
      mode_(block.mode()),
      variance_sampler_(prior.variance().parameter().prior()),
      block_(block),
      residual_(residual),
      normal_(0.0),
      rng_(rng)
{
    (void)std::get<bayes::SinglePerMarkerGaussianState>(block_.prior_state());
}

auto SinglePerMarkerGaussianStep::step() -> void
{
    auto& prior_state
        = std::get<bayes::SinglePerMarkerGaussianState>(block_.prior_state());
    auto& variance = prior_state.variance();
    auto& state = block_.state();
    const auto& xtx_diag = design_.xtx_diag(mode_);
    auto& coeffs = state.coeffs;

    normal_.reset();
    variance_sampler_.reset();
    {
        GeneticSweepAdjustmentGuard sweep{residual_, state};
        for (const Eigen::Index i : design_.valid_indices(mode_))
        {
            const double old_i = coeffs(i);
            GeneticResidualAdjustmentGuard guard{
                design_, mode_, i, state, residual_};
            const double rhs = design_.dot(mode_, i, residual_.y_adj)
                               + (xtx_diag(i) * old_i);
            normal_.set_prior_var(variance(i));
            coeffs(i) = normal_(
                NormalSampler<double>::Kernel{
                    .quadratic = xtx_diag(i),
                    .linear = rhs,
                    .scale = residual_.variance,
                },
                rng_);
            const double coeff = coeffs(i);
            variance(i) = variance_sampler_({1, coeff * coeff}, rng_);
        }
    }
    state.variance = detail::vecvar(state.u, detail::VarNormType::Population);
}

SingleSharedSpikeSlabStep::SingleSharedSpikeSlabStep(
    const bayes::GeneticDesign& design,
    const bayes::SingleSharedSpikeSlabGaussianPrior& prior,
    bayes::SingleGeneticBlockState& block,
    bayes::ResidualState& residual,
    std::mt19937_64& rng)
    : design_(design),
      mode_(block.mode()),
      variance_sampler_(prior.variance().parameter().prior()),
      proportion_sampler_(
          [&prior] -> std::optional<DirichletSampler<double>>
          {
              if (!prior.proportion().prior())
              {
                  return std::nullopt;
              }
              return DirichletSampler<double>{
                  prior.proportion().prior()->concentration()};
          }()),
      block_(block),
      residual_(residual),
      normal_(0.0),
      proportion_count_(
          Eigen::VectorXi::Zero(
              std::get<bayes::SingleSharedSpikeSlabGaussianState>(
                  block.prior_state())
                  .proportion()
                  .size())),
      rng_(rng)
{
}

auto SingleSharedSpikeSlabStep::step() -> void
{
    auto& prior_state = std::get<bayes::SingleSharedSpikeSlabGaussianState>(
        block_.prior_state());
    auto& variance = prior_state.variance();
    auto& assignment = prior_state.assignment();
    auto& proportion = prior_state.proportion();
    auto& state = block_.state();
    const auto& xtx_diag = design_.xtx_diag(mode_);
    auto& coeffs = state.coeffs;

    normal_.reset();
    normal_.set_prior_var(variance);
    uniform_.reset();
    variance_sampler_.reset();
    if (proportion_sampler_)
    {
        proportion_sampler_->reset();
    }
    logpi_ = proportion.array().log();
    proportion_count_.setZero();

    Eigen::Index variance_n = 0;
    double sum_squares = 0.0;
    {
        GeneticSweepAdjustmentGuard sweep{residual_, state};
        for (const Eigen::Index i : design_.valid_indices(mode_))
        {
            const double old_i = coeffs(i);
            const double rhs = design_.dot(mode_, i, residual_.y_adj)
                               + (xtx_diag(i) * old_i);
            const auto post = normal_.posterior_with_logL(
                NormalSampler<double>::Kernel{
                    .quadratic = xtx_diag(i),
                    .linear = rhs,
                    .scale = residual_.variance,
                });
            const double log_like_1_minus_0
                = post.log_likelihood_kernel + logpi_(1) - logpi_(0);
            const double prob_component_0
                = 1.0 / (1.0 + std::exp(log_like_1_minus_0));
            const int component = uniform_(rng_) < prob_component_0 ? 0 : 1;

            GeneticResidualAdjustmentGuard guard{
                design_, mode_, i, state, residual_};
            coeffs(i) = component == 0 ? 0.0 : normal_.draw(post.params, rng_);
            assignment(i) = component;

            ++proportion_count_(component);
            if (component != 0)
            {
                ++variance_n;
                const double coeff = coeffs(i);
                sum_squares += coeff * coeff;
            }
        }
    }

    variance = variance_sampler_({variance_n, sum_squares}, rng_);
    if (proportion_sampler_)
    {
        proportion = (*proportion_sampler_)(proportion_count_, rng_);
    }
    state.variance = detail::vecvar(state.u, detail::VarNormType::Population);
}

SinglePerMarkerSpikeSlabStep::SinglePerMarkerSpikeSlabStep(
    const bayes::GeneticDesign& design,
    const bayes::SinglePerMarkerSpikeSlabGaussianPrior& prior,
    bayes::SingleGeneticBlockState& block,
    bayes::ResidualState& residual,
    std::mt19937_64& rng)
    : design_(design),
      mode_(block.mode()),
      variance_sampler_(prior.variance().parameter().prior()),
      proportion_sampler_(
          [&prior] -> std::optional<DirichletSampler<double>>
          {
              if (!prior.proportion().prior())
              {
                  return std::nullopt;
              }
              return DirichletSampler<double>{
                  prior.proportion().prior()->concentration()};
          }()),
      block_(block),
      residual_(residual),
      normal_(0.0),
      proportion_count_(
          Eigen::VectorXi::Zero(
              std::get<bayes::SinglePerMarkerSpikeSlabGaussianState>(
                  block.prior_state())
                  .proportion()
                  .size())),
      rng_(rng)
{
}

auto SinglePerMarkerSpikeSlabStep::step() -> void
{
    auto& prior_state = std::get<bayes::SinglePerMarkerSpikeSlabGaussianState>(
        block_.prior_state());
    auto& variance = prior_state.variance();
    auto& assignment = prior_state.assignment();
    auto& proportion = prior_state.proportion();
    auto& state = block_.state();
    const auto& xtx_diag = design_.xtx_diag(mode_);
    auto& coeffs = state.coeffs;

    normal_.reset();
    uniform_.reset();
    variance_sampler_.reset();
    if (proportion_sampler_)
    {
        proportion_sampler_->reset();
    }
    logpi_ = proportion.array().log();
    proportion_count_.setZero();

    {
        GeneticSweepAdjustmentGuard sweep{residual_, state};
        for (const Eigen::Index i : design_.valid_indices(mode_))
        {
            const double old_i = coeffs(i);
            const double rhs = design_.dot(mode_, i, residual_.y_adj)
                               + (xtx_diag(i) * old_i);
            const auto post = normal_.set_prior_var(variance(i))
                                  .posterior_with_logL(
                                      NormalSampler<double>::Kernel{
                                          .quadratic = xtx_diag(i),
                                          .linear = rhs,
                                          .scale = residual_.variance,
                                      });
            const double log_like_1_minus_0
                = post.log_likelihood_kernel + logpi_(1) - logpi_(0);
            const double prob_component_0
                = 1.0 / (1.0 + std::exp(log_like_1_minus_0));
            const int component = uniform_(rng_) < prob_component_0 ? 0 : 1;

            GeneticResidualAdjustmentGuard guard{
                design_, mode_, i, state, residual_};
            coeffs(i) = component == 0 ? 0.0 : normal_.draw(post.params, rng_);
            assignment(i) = component;

            ++proportion_count_(component);
            if (component != 0)
            {
                const double coeff = coeffs(i);
                variance(i) = variance_sampler_({1, coeff * coeff}, rng_);
            }
        }
    }

    if (proportion_sampler_)
    {
        proportion = (*proportion_sampler_)(proportion_count_, rng_);
    }
    state.variance = detail::vecvar(state.u, detail::VarNormType::Population);
}

SingleScaledMixtureStep::SingleScaledMixtureStep(
    const bayes::GeneticDesign& design,
    const bayes::SingleScaledMixtureGaussianPrior& prior,
    bayes::SingleGeneticBlockState& block,
    bayes::ResidualState& residual,
    std::mt19937_64& rng)
    : design_(design),
      mode_(block.mode()),
      variance_sampler_(prior.variance().parameter().prior()),
      proportion_sampler_(prior.proportion().prior()->concentration()),
      multiplier_(prior.multiplier()),
      block_(block),
      residual_(residual),
      normal_(0.0),
      rng_(rng)
{
    assert(
        multiplier_.size() <= MAX_MIXTURE_COMPONENTS
        && "SingleScaledMixtureStep: mixture components exceed "
           "MAX_MIXTURE_COMPONENTS");
    marker_variances_.resize(multiplier_.size());
    logpi_.resize(multiplier_.size());
    const auto& prior_state = std::get<bayes::SingleScaledMixtureGaussianState>(
        block_.prior_state());
    proportion_count_ = Eigen::VectorXi::Zero(prior_state.proportion().size());
}

auto SingleScaledMixtureStep::step() -> void
{
    auto& prior_state = std::get<bayes::SingleScaledMixtureGaussianState>(
        block_.prior_state());
    auto& variance = prior_state.variance();
    auto& assignment = prior_state.assignment();
    auto& proportion = prior_state.proportion();
    auto& component = prior_state.component();
    auto& state = block_.state();
    const auto& xtx_diag = design_.xtx_diag(mode_);
    auto& coeffs = state.coeffs;

    normal_.reset();
    uniform_.reset();
    variance_sampler_.reset();
    proportion_sampler_.reset();
    logpi_ = proportion.array().log();
    marker_variances_ = variance * multiplier_.array();
    proportion_count_.setZero();

    Eigen::Index variance_n = 0;
    double sum_squares = 0.0;
    {
        GeneticSweepAdjustmentGuard sweep{residual_, state};
        for (const Eigen::Index i : design_.valid_indices(mode_))
        {
            const double old_i = coeffs(i);
            const double rhs = design_.dot(mode_, i, residual_.y_adj)
                               + (xtx_diag(i) * old_i);

            const Eigen::Index num_components = multiplier_.size();
            Eigen::Array<double, MAX_MIXTURE_COMPONENTS, 1> ll;
            ll(0) = logpi_(0);
            double max_ll = ll(0);
            for (Eigen::Index cls = 1; cls < num_components; ++cls)
            {
                scale_posts_[cls]
                    = normal_.set_prior_var(marker_variances_(cls))
                          .posterior_with_logL(
                              {.quadratic = xtx_diag(i),
                               .linear = rhs,
                               .scale = residual_.variance});
                ll(cls) = scale_posts_[cls].log_likelihood_kernel + logpi_(cls);
                max_ll = std::max(max_ll, ll(cls));
            }

            Eigen::Array<double, MAX_MIXTURE_COMPONENTS, 1> probs;
            double total = 0.0;
            for (Eigen::Index cls = 0; cls < num_components; ++cls)
            {
                probs(cls) = std::exp(ll(cls) - max_ll);
                total += probs(cls);
            }

            const double threshold = uniform_(rng_) * total;
            int component = static_cast<int>(num_components - 1);
            double cumsum = 0.0;
            for (Eigen::Index cls = 0; cls < num_components; ++cls)
            {
                cumsum += probs(cls);
                if (threshold < cumsum)
                {
                    component = static_cast<int>(cls);
                    break;
                }
            }

            GeneticResidualAdjustmentGuard residual_guard{
                design_, mode_, i, state, residual_};
            ComponentGebvAdjustmentGuard component_guard{
                design_, mode_, i, state, prior_state};
            coeffs(i) = 0.0;
            if (component > 0)
            {
                coeffs(i) = normal_.draw(scale_posts_[component].params, rng_);
            }
            assignment(i) = component;

            ++proportion_count_(component);
            if (component != 0)
            {
                ++variance_n;
                const double coeff = coeffs(i);
                sum_squares += (coeff * coeff) / multiplier_(component);
            }
        }
    }

    variance = variance_sampler_({variance_n, sum_squares}, rng_);
    proportion = proportion_sampler_(proportion_count_, rng_);
    state.variance = detail::vecvar(state.u, detail::VarNormType::Population);
    for (Eigen::Index k = 0;
         k < static_cast<Eigen::Index>(component.gebv.size());
         ++k)
    {
        component.gebv_var(k) = detail::vecvar(
            component.gebv[k], detail::VarNormType::Population);
    }
}

}  // namespace gelex
