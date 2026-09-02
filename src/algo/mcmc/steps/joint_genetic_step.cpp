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

#include "gelex/algo/mcmc/steps/joint_genetic_step.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <optional>
#include <random>
#include <utility>
#include <variant>

#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/half_normal_prior.h"
#include "gelex/bayes/genetic/half_normal_prior_state.h"
#include "gelex/bayes/genetic/legacy_genetic_prior.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/stats/beta_sampler.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/infra/stats/dirichlet_sampler.h"
#include "gelex/infra/stats/half_normal_sampler.h"
#include "gelex/infra/stats/normal_sampler.h"
#include "gelex/infra/stats/scaled_inv_chi2_sampler.h"
#include "gelex/types/genetic_mode.h"

#include "algo/mcmc/invariant.h"

namespace gelex
{

JointGaussianMixtureStep::JointGaussianMixtureStep(
    const bayes::GeneticDesign& design,
    const bayes::JointGaussianMixturePrior& prior,
    bayes::JointGeneticBlockState& block,
    bayes::ResidualState& residual,
    std::mt19937_64& rng)
    : design_(design),
      variance_samplers_(
          [&prior]
          {
              return std::array{
                  ScaledInvChi2Sampler<double>{
                      prior.variance(GeneticMode::A).parameter().prior()},
                  ScaledInvChi2Sampler<double>{
                      prior.variance(GeneticMode::D).parameter().prior()}};
          }()),
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
              std::get<bayes::JointGaussianMixtureState>(block.prior_state())
                  .proportion()
                  .size())),
      rng_(rng)
{
    std::ranges::set_intersection(
        design_.valid_indices(GeneticMode::A),
        design_.valid_indices(GeneticMode::D),
        std::back_inserter(valid_indices_));
}

auto JointGaussianMixtureStep::step() -> void
{
    auto& prior_state
        = std::get<bayes::JointGaussianMixtureState>(block_.prior_state());
    auto& additive_variance = prior_state.variance(GeneticMode::A);
    auto& dominance_variance = prior_state.variance(GeneticMode::D);
    auto& assignment = prior_state.assignment();
    auto& proportion = prior_state.proportion();
    auto& additive = block_.state(GeneticMode::A);
    auto& dominance = block_.state(GeneticMode::D);
    const auto& additive_xtx_diag = design_.xtx_diag(GeneticMode::A);
    const auto& dominance_xtx_diag = design_.xtx_diag(GeneticMode::D);
    auto& additive_coeffs = additive.coeffs;
    auto& dominance_coeffs = dominance.coeffs;

    normal_.reset();
    uniform_.reset();
    for (auto& sampler : variance_samplers_)
    {
        sampler.reset();
    }
    if (proportion_sampler_)
    {
        proportion_sampler_->reset();
    }
    logpi_ = proportion.array().log();
    proportion_count_.setZero();

    std::array<Eigen::Index, 2> variance_n{0, 0};
    std::array<double, 2> sum_squares{0.0, 0.0};
    for (const Eigen::Index i : valid_indices_)
    {
        const double old_additive_i = additive_coeffs(i);
        const double old_dominance_i = dominance_coeffs(i);
        const double additive_rhs
            = design_.dot(GeneticMode::A, i, residual_.y_adj)
              + (additive_xtx_diag(i) * old_additive_i);
        const double dominance_rhs
            = design_.dot(GeneticMode::D, i, residual_.y_adj)
              + (dominance_xtx_diag(i) * old_dominance_i);
        const auto additive_post
            = normal_.set_prior_var(additive_variance)
                  .posterior_with_logL(
                      NormalSampler<double>::Kernel{
                          .quadratic = additive_xtx_diag(i),
                          .linear = additive_rhs,
                          .scale = residual_.variance,
                      });
        const auto dominance_post
            = normal_.set_prior_var(dominance_variance)
                  .posterior_with_logL(
                      NormalSampler<double>::Kernel{
                          .quadratic = dominance_xtx_diag(i),
                          .linear = dominance_rhs,
                          .scale = residual_.variance,
                      });

        Eigen::Array<double, 4, 1> log_likelihoods;
        log_likelihoods(0) = logpi_(0);
        log_likelihoods(1) = additive_post.log_likelihood_kernel + logpi_(1);
        log_likelihoods(2) = dominance_post.log_likelihood_kernel + logpi_(2);
        log_likelihoods(3) = additive_post.log_likelihood_kernel
                             + dominance_post.log_likelihood_kernel + logpi_(3);
        const double max_log_likelihood = log_likelihoods.maxCoeff();
        const auto probabilities = (log_likelihoods - max_log_likelihood).exp();
        const double total = probabilities.sum();
        const double threshold = uniform_(rng_) * total;

        int component = 3;
        double cumsum = 0.0;
        for (Eigen::Index cls = 0; cls < probabilities.size(); ++cls)
        {
            cumsum += probabilities(cls);
            if (threshold < cumsum)
            {
                component = static_cast<int>(cls);
                break;
            }
        }

        JointGeneticAdjustmentGuard guard{design_, i, block_, residual_};
        additive_coeffs(i) = (component == 1 || component == 3)
                                 ? normal_.draw(additive_post.params, rng_)
                                 : 0.0;
        dominance_coeffs(i) = (component == 2 || component == 3)
                                  ? normal_.draw(dominance_post.params, rng_)
                                  : 0.0;
        assignment(i) = component;

        ++proportion_count_(component);
        if (component == 1 || component == 3)
        {
            ++variance_n[std::to_underlying(GeneticMode::A)];
            const double coeff = additive_coeffs(i);
            sum_squares[std::to_underlying(GeneticMode::A)] += coeff * coeff;
        }
        if (component == 2 || component == 3)
        {
            ++variance_n[std::to_underlying(GeneticMode::D)];
            const double coeff = dominance_coeffs(i);
            sum_squares[std::to_underlying(GeneticMode::D)] += coeff * coeff;
        }
    }

    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (auto mode : modes)
    {
        const auto index = std::to_underlying(mode);
        prior_state.variance(mode) = variance_samplers_[index](
            {variance_n[index], sum_squares[index]}, rng_);
    }
    if (proportion_sampler_)
    {
        proportion = (*proportion_sampler_)(proportion_count_, rng_);
    }
    additive.variance
        = detail::vecvar(additive.u, detail::VarNormType::Population);
    dominance.variance
        = detail::vecvar(dominance.u, detail::VarNormType::Population);
}

JointHalfNormalMixtureStep::JointHalfNormalMixtureStep(
    const bayes::GeneticDesign& design,
    const bayes::JointHalfNormalMixturePrior& prior,
    bayes::JointGeneticBlockState& block,
    bayes::ResidualState& residual,
    std::mt19937_64& rng)
    : design_(design),
      variance_samplers_(
          [&prior]
          {
              return std::array{
                  ScaledInvChi2Sampler<double>{
                      prior.variance(GeneticMode::A).parameter().prior()},
                  ScaledInvChi2Sampler<double>{
                      prior.variance(GeneticMode::D).parameter().prior()}};
          }()),
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
      half_normal_(
          std::get<bayes::JointHalfNormalMixtureState>(block.prior_state())
              .variance(GeneticMode::D)),
      sign_sampler_(
          prior.dominance_positive_probability().prior().alpha(),
          prior.dominance_positive_probability().prior().beta()),
      proportion_count_(
          Eigen::VectorXi::Zero(
              std::get<bayes::JointHalfNormalMixtureState>(block.prior_state())
                  .proportion()
                  .size())),
      rng_(rng)
{
    std::ranges::set_intersection(
        design_.valid_indices(GeneticMode::A),
        design_.valid_indices(GeneticMode::D),
        std::back_inserter(valid_indices_));
}

auto JointHalfNormalMixtureStep::step() -> void
{
    auto& prior_state
        = std::get<bayes::JointHalfNormalMixtureState>(block_.prior_state());
    auto& additive_variance = prior_state.variance(GeneticMode::A);
    auto& dominance_variance = prior_state.variance(GeneticMode::D);
    auto& assignment = prior_state.assignment();
    auto& proportion = prior_state.proportion();
    auto& dominance_sign = prior_state.dominance_sign();
    auto& additive = block_.state(GeneticMode::A);
    auto& dominance = block_.state(GeneticMode::D);
    const auto& additive_xtx_diag = design_.xtx_diag(GeneticMode::A);
    const auto& dominance_xtx_diag = design_.xtx_diag(GeneticMode::D);
    auto& additive_coeffs = additive.coeffs;
    auto& dominance_coeffs = dominance.coeffs;

    normal_.reset();
    half_normal_.reset();
    uniform_.reset();
    sign_sampler_.reset();
    for (auto& sampler : variance_samplers_)
    {
        sampler.reset();
    }
    if (proportion_sampler_)
    {
        proportion_sampler_->reset();
    }
    logpi_ = proportion.array().log();
    proportion_count_.setZero();

    std::array<Eigen::Index, 2> variance_n{0, 0};
    std::array<double, 2> sum_squares{0.0, 0.0};
    BetaSampler<double>::Likelihood sign_likelihood{};
    for (const Eigen::Index i : valid_indices_)
    {
        const double old_additive_i = additive_coeffs(i);
        const double old_dominance_i = dominance_coeffs(i);
        const double additive_rhs
            = design_.dot(GeneticMode::A, i, residual_.y_adj)
              + (additive_xtx_diag(i) * old_additive_i);
        const double dominance_rhs
            = design_.dot(GeneticMode::D, i, residual_.y_adj)
              + (dominance_xtx_diag(i) * old_dominance_i);
        const auto additive_post
            = normal_.set_prior_var(additive_variance)
                  .posterior_with_logL(
                      NormalSampler<double>::Kernel{
                          .quadratic = additive_xtx_diag(i),
                          .linear = additive_rhs,
                          .scale = residual_.variance,
                      });
        const auto dominance_positive_post
            = half_normal_.set_prior_var(dominance_variance)
                  .posterior_with_logL(
                      HalfNormalSampler<double>::Kernel{
                          .quadratic = dominance_xtx_diag(i),
                          .linear = dominance_rhs,
                          .scale = residual_.variance,
                      },
                      static_cast<std::int8_t>(1));
        const auto dominance_negative_post = half_normal_.posterior_with_logL(
            HalfNormalSampler<double>::Kernel{
                .quadratic = dominance_xtx_diag(i),
                .linear = dominance_rhs,
                .scale = residual_.variance,
            },
            static_cast<std::int8_t>(-1));

        const double positive_sign_log_weight
            = std::log(dominance_sign.positive_probability)
              + dominance_positive_post.log_marginal_kernel;
        const double negative_sign_log_weight
            = std::log(1.0 - dominance_sign.positive_probability)
              + dominance_negative_post.log_marginal_kernel;
        const double max_sign_log_weight
            = positive_sign_log_weight > negative_sign_log_weight
                  ? positive_sign_log_weight
                  : negative_sign_log_weight;
        const double dominance_log_likelihood
            = max_sign_log_weight
              + std::log(
                  std::exp(positive_sign_log_weight - max_sign_log_weight)
                  + std::exp(negative_sign_log_weight - max_sign_log_weight));

        Eigen::Array<double, 4, 1> log_likelihoods;
        log_likelihoods(0) = logpi_(0);
        log_likelihoods(1) = additive_post.log_likelihood_kernel + logpi_(1);
        log_likelihoods(2) = dominance_log_likelihood + logpi_(2);
        log_likelihoods(3) = additive_post.log_likelihood_kernel
                             + dominance_log_likelihood + logpi_(3);
        const double max_log_likelihood = log_likelihoods.maxCoeff();
        const auto probabilities = (log_likelihoods - max_log_likelihood).exp();
        const double total = probabilities.sum();
        const double threshold = uniform_(rng_) * total;

        int component = 3;
        double cumsum = 0.0;
        for (Eigen::Index cls = 0; cls < probabilities.size(); ++cls)
        {
            cumsum += probabilities(cls);
            if (threshold < cumsum)
            {
                component = static_cast<int>(cls);
                break;
            }
        }

        const bool dominance_active = component == 2 || component == 3;
        const auto* dominance_post = &dominance_negative_post;
        if (dominance_active)
        {
            const double sign_total
                = std::exp(positive_sign_log_weight - max_sign_log_weight)
                  + std::exp(negative_sign_log_weight - max_sign_log_weight);
            const double sign_threshold = uniform_(rng_) * sign_total;
            const double positive_sign_probability
                = std::exp(positive_sign_log_weight - max_sign_log_weight);
            if (sign_threshold < positive_sign_probability)
            {
                dominance_sign.sign(i) = 1;
                dominance_post = &dominance_positive_post;
                ++sign_likelihood.n_success;
            }
            else
            {
                dominance_sign.sign(i) = 0;
                ++sign_likelihood.n_fail;
            }
        }

        JointGeneticAdjustmentGuard guard{design_, i, block_, residual_};
        additive_coeffs(i) = (component == 1 || component == 3)
                                 ? normal_.draw(additive_post.params, rng_)
                                 : 0.0;
        dominance_coeffs(i)
            = dominance_active ? half_normal_.draw(*dominance_post, rng_) : 0.0;
        assignment(i) = component;

        ++proportion_count_(component);
        if (component == 1 || component == 3)
        {
            ++variance_n[std::to_underlying(GeneticMode::A)];
            const double coeff = additive_coeffs(i);
            sum_squares[std::to_underlying(GeneticMode::A)] += coeff * coeff;
        }
        if (dominance_active)
        {
            ++variance_n[std::to_underlying(GeneticMode::D)];
            const double coeff = dominance_coeffs(i);
            sum_squares[std::to_underlying(GeneticMode::D)] += coeff * coeff;
        }
    }

    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    for (auto mode : modes)
    {
        const auto index = std::to_underlying(mode);
        prior_state.variance(mode) = variance_samplers_[index](
            {variance_n[index], sum_squares[index]}, rng_);
    }
    if (proportion_sampler_)
    {
        proportion = (*proportion_sampler_)(proportion_count_, rng_);
    }
    dominance_sign.positive_probability = sign_sampler_(sign_likelihood, rng_);
    additive.variance
        = detail::vecvar(additive.u, detail::VarNormType::Population);
    dominance.variance
        = detail::vecvar(dominance.u, detail::VarNormType::Population);
}

}  // namespace gelex
