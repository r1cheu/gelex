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

#ifndef GELEX_ALGO_INFER_MCMC_TRANSITIONS_GENETIC_COEFFICIENT_H_
#define GELEX_ALGO_INFER_MCMC_TRANSITIONS_GENETIC_COEFFICIENT_H_

#include <cmath>
#include <random>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/invariant.h"
#include "gelex/algo/infer/mcmc/kernels/detail/mixture_op.h"
#include "gelex/exception.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class SingleSharedGaussianCoefficientTransition
{
   public:
    SingleSharedGaussianCoefficientTransition(
        const bayes::SingleGeneticPrior& /*prior*/,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual)
        : state_(block.state()),
          residual_(residual),
          variance_(block.prior_state()
                        .require<bayes::SingleSharedVarianceStateCap>()
                        .variance()),
          normal_(variance_)
    {
    }

    auto prepare() -> void { normal_.reset(); }

    auto update(
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

    auto finish() -> void { state_.variance = stats::detail::var(state_.u)(0); }

   private:
    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;
    double& variance_;
    stats::NormalSampler<double> normal_;
};

class SinglePerMarkerGaussianCoefficientTransition
{
   public:
    SinglePerMarkerGaussianCoefficientTransition(
        const bayes::SingleGeneticPrior& /*prior*/,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual)
        : state_(block.state()),
          residual_(residual),
          variance_(block.prior_state()
                        .require<bayes::SinglePerMarkerVarianceStateCap>()
                        .variance()),
          normal_(0.0)
    {
    }

    auto prepare() -> void { normal_.reset(); }

    auto update(
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

    auto finish() -> void { state_.variance = stats::detail::var(state_.u)(0); }

   private:
    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;
    Eigen::VectorXd& variance_;
    stats::NormalSampler<double> normal_;
};

class SingleSharedSpikeSlabCoefficientTransition
{
   public:
    SingleSharedSpikeSlabCoefficientTransition(
        const bayes::SingleGeneticPrior& /*prior*/,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual)
        : state_(block.state()),
          residual_(residual),
          variance_(block.prior_state()
                        .require<bayes::SingleSharedVarianceStateCap>()
                        .variance()),
          proportion_(block.prior_state()
                          .require<bayes::SingleProportionStateCap>()
                          .proportion()),
          normal_(variance_)
    {
    }

    auto prepare() -> void
    {
        normal_.reset();
        uniform_.reset();
        logpi_ = proportion_.value.array().log();
    }

    auto update(
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
        const double prob_component_0
            = 1.0 / (1.0 + std::exp(log_like_1_minus_0));
        const int component = uniform_(rng) < prob_component_0 ? 0 : 1;

        ProportionAssignmentGuard assignment_guard{proportion_, marker_index};
        GeneticAdjustmentGuard guard{
            column, coeffs(marker_index), residual_, state_};
        coeffs(marker_index)
            = component == 0 ? 0.0 : normal_.draw(post.params, rng);
        proportion_.assignment(marker_index) = component;
    }

    auto finish() -> void { state_.variance = stats::detail::var(state_.u)(0); }

   private:
    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;
    double& variance_;
    bayes::ProportionState& proportion_;
    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
};

class SinglePerMarkerSpikeSlabCoefficientTransition
{
   public:
    SinglePerMarkerSpikeSlabCoefficientTransition(
        const bayes::SingleGeneticPrior& /*prior*/,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual)
        : state_(block.state()),
          residual_(residual),
          variance_(block.prior_state()
                        .require<bayes::SinglePerMarkerVarianceStateCap>()
                        .variance()),
          proportion_(block.prior_state()
                          .require<bayes::SingleProportionStateCap>()
                          .proportion()),
          normal_(0.0)
    {
    }

    auto prepare() -> void
    {
        normal_.reset();
        uniform_.reset();
        logpi_ = proportion_.value.array().log();
    }

    auto update(
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
        const double prob_component_0
            = 1.0 / (1.0 + std::exp(log_like_1_minus_0));
        const int component = uniform_(rng) < prob_component_0 ? 0 : 1;

        ProportionAssignmentGuard assignment_guard{proportion_, marker_index};
        GeneticAdjustmentGuard guard{
            column, coeffs(marker_index), residual_, state_};
        coeffs(marker_index)
            = component == 0 ? 0.0 : normal_.draw(post.params, rng);
        proportion_.assignment(marker_index) = component;
    }

    auto finish() -> void { state_.variance = stats::detail::var(state_.u)(0); }

   private:
    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;
    Eigen::VectorXd& variance_;
    bayes::ProportionState& proportion_;
    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    Eigen::VectorXd logpi_;
};

class SingleScaledMixtureCoefficientTransition
{
   public:
    SingleScaledMixtureCoefficientTransition(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        bayes::ResidualState& residual)
        : state_(block.state()),
          residual_(residual),
          variance_(block.prior_state()
                        .require<bayes::SingleSharedVarianceStateCap>()
                        .variance()),
          proportion_(block.prior_state()
                          .require<bayes::SingleProportionStateCap>()
                          .proportion()),
          component_(block.prior_state()
                         .require<bayes::SingleComponentStateCap>()
                         .component()),
          normal_(0.0)
    {
        const auto* multiplier_cap = prior.query<bayes::SingleMultiplierCap>();
        if (multiplier_cap == nullptr)
        {
            throw GelexException(
                "SingleScaledMixtureCoefficientTransition: prior lacks "
                "multiplier capability");
        }
        multiplier_ = multiplier_cap->multiplier();
        marker_variances_.resize(multiplier_.size());
        logpi_.resize(multiplier_.size());
    }

    auto prepare() -> void
    {
        normal_.reset();
        uniform_.reset();
        logpi_ = proportion_.value.array().log();
        marker_variances_ = variance_ * multiplier_.array();
    }

    auto update(
        Eigen::Index marker_index,
        Eigen::Ref<const Eigen::VectorXd> column,
        double xtx_diag_i,
        std::mt19937_64& rng) -> void
    {
        auto& coeffs = state_.coeffs;
        const double old_i = coeffs(marker_index);
        const double rhs = column.dot(residual_.y_adj) + (xtx_diag_i * old_i);
        const Eigen::Index num_components = multiplier_.size();
        scale_pp_.log_likelihoods(0) = 0.0;
        for (Eigen::Index cls = 1; cls < num_components; ++cls)
        {
            const auto post = normal_.set_prior_var(marker_variances_(cls))
                                  .posterior_with_logL(
                                      {.quadratic = xtx_diag_i,
                                       .linear = rhs,
                                       .scale = residual_.variance});
            scale_pp_.means(cls) = post.params.mean;
            scale_pp_.vars(cls) = post.params.var;
            scale_pp_.log_likelihoods(cls) = post.log_likelihood_kernel;
        }

        Eigen::Array<double, detail::kMaxMixtureComponents, 1> ll;
        Eigen::Array<double, detail::kMaxMixtureComponents, 1> probs;
        ll.head(num_components) = scale_pp_.log_likelihoods.head(num_components)
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
                {.mean = scale_pp_.means(component),
                 .var = scale_pp_.vars(component)},
                rng);
        }
        proportion_.assignment(marker_index) = component;
    }

    auto finish() -> void
    {
        state_.variance = stats::detail::var(state_.u)(0);
        detail::compute_component_variances(component_);
    }

   private:
    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;
    double& variance_;
    bayes::ProportionState& proportion_;
    bayes::ComponentState& component_;
    Eigen::VectorXd multiplier_;
    Eigen::VectorXd marker_variances_;
    Eigen::VectorXd logpi_;
    stats::NormalSampler<double> normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    detail::MixtureNormalPosteriors scale_pp_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_TRANSITIONS_GENETIC_COEFFICIENT_H_
