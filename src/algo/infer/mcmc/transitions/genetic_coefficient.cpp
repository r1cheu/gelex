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

#include "gelex/algo/infer/mcmc/transitions/genetic_coefficient.h"

#include <cmath>

#include "gelex/algo/infer/mcmc/invariant.h"
#include "gelex/infra/stats/detail/var.h"

namespace gelex::mcmc
{

SingleSharedGaussianTransition::SingleSharedGaussianTransition(
    const bayes::SingleGeneticPrior& /*prior*/,
    bayes::SingleGeneticBlockState& block,
    bayes::ResidualState& residual)
    : state_(block.state()),
      residual_(residual),
      variance_(block.prior_state()
                    .get<bayes::SingleSharedVarianceStateCap>()
                    .variance()),
      normal_(variance_)
{
}

auto SingleSharedGaussianTransition::prepare() -> void
{
    normal_.reset();
}

auto SingleSharedGaussianTransition::update(
    Eigen::Index marker_index,
    Eigen::Ref<const Eigen::VectorXd> column,
    double xtx_diag_i,
    std::mt19937_64& rng) -> void
{
    auto& coeffs = state_.coeffs;
    const double old_i = coeffs(marker_index);
    GeneticAdjustmentGuard guard{
        column, coeffs(marker_index), residual_, state_};
    const double rhs = column.dot(residual_.y_adj) + (xtx_diag_i * old_i);
    normal_.set_prior_var(variance_);
    coeffs(marker_index) = normal_(
        stats::NormalSampler<double>::Kernel{
            .quadratic = xtx_diag_i,
            .linear = rhs,
            .scale = residual_.variance,
        },
        rng);
}

auto SingleSharedGaussianTransition::finish() -> void
{
    state_.variance = stats::detail::var(state_.u)(0);
}

SinglePerMarkerGaussianTransition::SinglePerMarkerGaussianTransition(
    const bayes::SingleGeneticPrior& /*prior*/,
    bayes::SingleGeneticBlockState& block,
    bayes::ResidualState& residual)
    : state_(block.state()),
      residual_(residual),
      variance_(block.prior_state()
                    .get<bayes::SinglePerMarkerVarianceStateCap>()
                    .variance()),
      normal_(0.0)
{
}

auto SinglePerMarkerGaussianTransition::prepare() -> void
{
    normal_.reset();
}

auto SinglePerMarkerGaussianTransition::update(
    Eigen::Index marker_index,
    Eigen::Ref<const Eigen::VectorXd> column,
    double xtx_diag_i,
    std::mt19937_64& rng) -> void
{
    auto& coeffs = state_.coeffs;
    const double old_i = coeffs(marker_index);
    GeneticAdjustmentGuard guard{
        column, coeffs(marker_index), residual_, state_};
    const double rhs = column.dot(residual_.y_adj) + (xtx_diag_i * old_i);
    normal_.set_prior_var(variance_(marker_index));
    coeffs(marker_index) = normal_(
        stats::NormalSampler<double>::Kernel{
            .quadratic = xtx_diag_i,
            .linear = rhs,
            .scale = residual_.variance,
        },
        rng);
}

auto SinglePerMarkerGaussianTransition::finish() -> void
{
    state_.variance = stats::detail::var(state_.u)(0);
}

SingleSharedSpikeSlabTransition::SingleSharedSpikeSlabTransition(
    const bayes::SingleGeneticPrior& /*prior*/,
    bayes::SingleGeneticBlockState& block,
    bayes::ResidualState& residual)
    : state_(block.state()),
      residual_(residual),
      variance_(block.prior_state()
                    .get<bayes::SingleSharedVarianceStateCap>()
                    .variance()),
      proportion_(block.prior_state()
                      .get<bayes::SingleProportionStateCap>()
                      .proportion()),
      normal_(variance_)
{
}

auto SingleSharedSpikeSlabTransition::prepare() -> void
{
    normal_.reset();
    uniform_.reset();
    logpi_ = proportion_.value.array().log();
}

auto SingleSharedSpikeSlabTransition::update(
    Eigen::Index marker_index,
    Eigen::Ref<const Eigen::VectorXd> column,
    double xtx_diag_i,
    std::mt19937_64& rng) -> void
{
    auto& coeffs = state_.coeffs;
    const double old_i = coeffs(marker_index);
    const double rhs = column.dot(residual_.y_adj) + (xtx_diag_i * old_i);
    const auto post = normal_.set_prior_var(variance_).posterior_with_logL(
        stats::NormalSampler<double>::Kernel{
            .quadratic = xtx_diag_i,
            .linear = rhs,
            .scale = residual_.variance,
        });
    const double log_like_1_minus_0
        = post.log_likelihood_kernel + logpi_(1) - logpi_(0);
    const double prob_component_0 = 1.0 / (1.0 + std::exp(log_like_1_minus_0));
    const int component = uniform_(rng) < prob_component_0 ? 0 : 1;

    ProportionAssignmentGuard assignment_guard{proportion_, marker_index};
    GeneticAdjustmentGuard guard{
        column, coeffs(marker_index), residual_, state_};
    coeffs(marker_index)
        = component == 0 ? 0.0 : normal_.draw(post.params, rng);
    proportion_.assignment(marker_index) = component;
}

auto SingleSharedSpikeSlabTransition::finish() -> void
{
    state_.variance = stats::detail::var(state_.u)(0);
}

SinglePerMarkerSpikeSlabTransition::SinglePerMarkerSpikeSlabTransition(
    const bayes::SingleGeneticPrior& /*prior*/,
    bayes::SingleGeneticBlockState& block,
    bayes::ResidualState& residual)
    : state_(block.state()),
      residual_(residual),
      variance_(block.prior_state()
                    .get<bayes::SinglePerMarkerVarianceStateCap>()
                    .variance()),
      proportion_(block.prior_state()
                      .get<bayes::SingleProportionStateCap>()
                      .proportion()),
      normal_(0.0)
{
}

auto SinglePerMarkerSpikeSlabTransition::prepare() -> void
{
    normal_.reset();
    uniform_.reset();
    logpi_ = proportion_.value.array().log();
}

auto SinglePerMarkerSpikeSlabTransition::update(
    Eigen::Index marker_index,
    Eigen::Ref<const Eigen::VectorXd> column,
    double xtx_diag_i,
    std::mt19937_64& rng) -> void
{
    auto& coeffs = state_.coeffs;
    const double old_i = coeffs(marker_index);
    const double rhs = column.dot(residual_.y_adj) + (xtx_diag_i * old_i);
    const auto post = normal_.set_prior_var(variance_(marker_index))
                          .posterior_with_logL(
                              stats::NormalSampler<double>::Kernel{
                                  .quadratic = xtx_diag_i,
                                  .linear = rhs,
                                  .scale = residual_.variance,
                              });
    const double log_like_1_minus_0
        = post.log_likelihood_kernel + logpi_(1) - logpi_(0);
    const double prob_component_0 = 1.0 / (1.0 + std::exp(log_like_1_minus_0));
    const int component = uniform_(rng) < prob_component_0 ? 0 : 1;

    ProportionAssignmentGuard assignment_guard{proportion_, marker_index};
    GeneticAdjustmentGuard guard{
        column, coeffs(marker_index), residual_, state_};
    coeffs(marker_index)
        = component == 0 ? 0.0 : normal_.draw(post.params, rng);
    proportion_.assignment(marker_index) = component;
}

auto SinglePerMarkerSpikeSlabTransition::finish() -> void
{
    state_.variance = stats::detail::var(state_.u)(0);
}

SingleScaledMixtureTransition::SingleScaledMixtureTransition(
    const bayes::SingleGeneticPrior& prior,
    bayes::SingleGeneticBlockState& block,
    bayes::ResidualState& residual)
    : state_(block.state()),
      residual_(residual),
      variance_(block.prior_state()
                    .get<bayes::SingleSharedVarianceStateCap>()
                    .variance()),
      proportion_(block.prior_state()
                      .get<bayes::SingleProportionStateCap>()
                      .proportion()),
      component_(block.prior_state()
                     .get<bayes::SingleComponentStateCap>()
                     .component()),
      normal_(0.0)
{
    multiplier_ = prior.get<bayes::SingleMultiplierCap>().multiplier();
    marker_variances_.resize(multiplier_.size());
    logpi_.resize(multiplier_.size());
}

auto SingleScaledMixtureTransition::prepare() -> void
{
    normal_.reset();
    uniform_.reset();
    logpi_ = proportion_.value.array().log();
    marker_variances_ = variance_ * multiplier_.array();
}

auto SingleScaledMixtureTransition::update(
    Eigen::Index marker_index,
    Eigen::Ref<const Eigen::VectorXd> column,
    double xtx_diag_i,
    std::mt19937_64& rng) -> void
{
    auto& coeffs = state_.coeffs;
    const double old_i = coeffs(marker_index);
    const double rhs = column.dot(residual_.y_adj) + (xtx_diag_i * old_i);
    const Eigen::Index num_components = multiplier_.size();
    scale_log_likelihoods_(0) = 0.0;
    for (Eigen::Index cls = 1; cls < num_components; ++cls)
    {
        const auto post = normal_.set_prior_var(marker_variances_(cls))
                              .posterior_with_logL(
                                  {.quadratic = xtx_diag_i,
                                   .linear = rhs,
                                   .scale = residual_.variance});
        scale_means_(cls) = post.params.mean;
        scale_vars_(cls) = post.params.var;
        scale_log_likelihoods_(cls) = post.log_likelihood_kernel;
    }

    Eigen::Array<double, kMaxMixtureComponents, 1> ll;
    Eigen::Array<double, kMaxMixtureComponents, 1> probs;
    ll.head(num_components) = scale_log_likelihoods_.head(num_components)
                              + logpi_.head(num_components).array();
    const double max_ll = ll.head(num_components).maxCoeff();
    probs.head(num_components) = (ll.head(num_components) - max_ll).exp();
    const double total = probs.head(num_components).sum();

    const double threshold = uniform_(rng) * total;
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

    ProportionAssignmentGuard assignment_guard{proportion_, marker_index};
    GeneticMixtureAdjustmentGuard guard{
        column,
        coeffs(marker_index),
        residual_,
        state_,
        component_,
        proportion_.assignment,
        marker_index};
    coeffs(marker_index) = 0.0;
    if (component > 0)
    {
        coeffs(marker_index) = normal_.draw(
            {.mean = scale_means_(component), .var = scale_vars_(component)},
            rng);
    }
    proportion_.assignment(marker_index) = component;
}

auto SingleScaledMixtureTransition::finish() -> void
{
    state_.variance = stats::detail::var(state_.u)(0);
    for (Eigen::Index k = 0;
         k < static_cast<Eigen::Index>(component_.gebv.size());
         ++k)
    {
        component_.gebv_var(k) = stats::detail::var(component_.gebv[k])(0);
    }
}

JointGaussianMixtureTransition::JointGaussianMixtureTransition(
    const bayes::JointGeneticPrior& /*prior*/,
    bayes::JointGeneticBlockState& block,
    bayes::ResidualState& residual)
    : additive_(block.state(GeneticMode::A)),
      dominance_(block.state(GeneticMode::D)),
      residual_(residual),
      variance_(block.prior_state().get<bayes::JointSharedVarianceStateCap>()),
      proportion_(block.prior_state()
                      .get<bayes::JointProportionStateCap>()
                      .proportion()),
      normal_(0.0)
{
}

auto JointGaussianMixtureTransition::prepare() -> void
{
    normal_.reset();
    uniform_.reset();
    logpi_ = proportion_.value.array().log();
}

auto JointGaussianMixtureTransition::update(
    Eigen::Index marker_index,
    Eigen::Ref<const Eigen::VectorXd> additive_column,
    double additive_xtx_diag_i,
    Eigen::Ref<const Eigen::VectorXd> dominance_column,
    double dominance_xtx_diag_i,
    std::mt19937_64& rng) -> void
{
    auto& additive_coeffs = additive_.coeffs;
    auto& dominance_coeffs = dominance_.coeffs;
    const double old_additive_i = additive_coeffs(marker_index);
    const double old_dominance_i = dominance_coeffs(marker_index);
    const double additive_rhs = additive_column.dot(residual_.y_adj)
                                + (additive_xtx_diag_i * old_additive_i);
    const double dominance_rhs = dominance_column.dot(residual_.y_adj)
                                 + (dominance_xtx_diag_i * old_dominance_i);
    const auto additive_post
        = normal_.set_prior_var(variance_.variance(GeneticMode::A))
              .posterior_with_logL(
                  stats::NormalSampler<double>::Kernel{
                      .quadratic = additive_xtx_diag_i,
                      .linear = additive_rhs,
                      .scale = residual_.variance,
                  });
    const auto dominance_post
        = normal_.set_prior_var(variance_.variance(GeneticMode::D))
              .posterior_with_logL(
                  stats::NormalSampler<double>::Kernel{
                      .quadratic = dominance_xtx_diag_i,
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
    const double threshold = uniform_(rng) * total;

    int component = 3;
    double cumsum = 0.0;
    for (Eigen::Index i = 0; i < probabilities.size(); ++i)
    {
        cumsum += probabilities(i);
        if (threshold < cumsum)
        {
            component = static_cast<int>(i);
            break;
        }
    }

    ProportionAssignmentGuard assignment_guard{proportion_, marker_index};
    JointGeneticAdjustmentGuard guard{
        additive_column,
        dominance_column,
        additive_coeffs(marker_index),
        dominance_coeffs(marker_index),
        residual_,
        additive_,
        dominance_};
    additive_coeffs(marker_index)
        = (component == 1 || component == 3)
              ? normal_.draw(additive_post.params, rng)
              : 0.0;
    dominance_coeffs(marker_index)
        = (component == 2 || component == 3)
              ? normal_.draw(dominance_post.params, rng)
              : 0.0;
    proportion_.assignment(marker_index) = component;
}

auto JointGaussianMixtureTransition::finish() -> void
{
    additive_.variance = stats::detail::var(additive_.u)(0);
    dominance_.variance = stats::detail::var(dominance_.u)(0);
}

}  // namespace gelex::mcmc
