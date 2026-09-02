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

#ifndef GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_KERNEL_H_
#define GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_KERNEL_H_

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>

#include "gelex/bayes/design.h"
#include "gelex/bayes/detail/fitted_value_update.h"
#include "gelex/bayes/genetic/joint_topology.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/probability_updater.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/stats/half_normal_sampler.h"
#include "gelex/infra/stats/normal_sampler.h"
#include "gelex/infra/stats/scaled_inv_chi2_sampler.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

template <
    UpdatePolicy ProbabilitiesUpdate,
    UpdatePolicy PositiveProbabilityUpdate>
class JointSpikeSlabKernel
{
    using ModePrior = GaussianPrior<VarianceLayout::Pooled>;
    using JointPrior
        = JointSpikeSlabPrior<ProbabilitiesUpdate, PositiveProbabilityUpdate>;
    using GeneticPrior = JointTopology<ModePrior, JointPrior>;
    using Layout = typename JointPrior::component_layout;
    using ModeState = GeneticModeState<
        GaussianState<VarianceLayout::Pooled>,
        typename JointPrior::component_layout>;
    using GeneticState = JointTopology<ModeState, JointSpikeSlabState>;
    using DominancePosterior = HalfNormalSampler<double>::Posterior;

    static constexpr std::size_t class_count = JointPrior::class_count;

    struct MarkerLikelihood
    {
        double additive_quadratic{};
        double additive_linear{};
        double dominance_quadratic{};
        double dominance_linear{};
        double residual_variance{};
    };

    struct MarkerSample
    {
        std::size_t class_index{};
        double additive{};
        double dominance{};
        bool positive_dominance{};
    };

   public:
    explicit JointSpikeSlabKernel(const GeneticPrior& prior)
        : additive_variance_sampler_{prior.mode_values()
                                         .template get<GeneticMode::A>()
                                         .variance.prior()},
          dominance_variance_sampler_{prior.mode_values()
                                          .template get<GeneticMode::D>()
                                          .variance.prior()},
          probabilities_updater_{prior.joint().probabilities},
          positive_probability_updater_{prior.joint().positive_probability},
          half_normal_{prior.mode_values()
                           .template get<GeneticMode::D>()
                           .variance.initial_value()}
    {
    }

    auto step(
        const bayes::GeneticDesign& design,
        GeneticState& state,
        ResidualState& residual,
        std::mt19937_64& rng) -> void
    {
        const auto& additive_projection = design.projection(GeneticMode::A);
        const auto& dominance_projection = design.projection(GeneticMode::D);
        const auto valid_indices = design.common_valid_indices();
        auto& additive = state.mode_values().template get<GeneticMode::A>();
        auto& dominance = state.mode_values().template get<GeneticMode::D>();
        auto& joint = state.joint();
        auto& additive_variance = additive.family_state.variance;
        auto& dominance_variance = dominance.family_state.variance;
        const auto& additive_xtx_diag = additive_projection.xtx_diag();
        const auto& dominance_xtx_diag = dominance_projection.xtx_diag();

        reset_samplers();
        if constexpr (ProbabilitiesUpdate == UpdatePolicy::Sampled)
        {
            probabilities_updater_.reset();
        }
        if constexpr (PositiveProbabilityUpdate == UpdatePolicy::Sampled)
        {
            positive_probability_updater_.reset();
        }

        std::array<double, class_count> log_probabilities{};
        for (std::size_t class_index = 0; class_index < class_count;
             ++class_index)
        {
            log_probabilities[class_index]
                = std::log(joint.probabilities[class_index]);
        }
        normal_.set_prior_var(additive_variance);
        half_normal_.set_prior_var(dominance_variance);

        Eigen::Index additive_active_count = 0;
        Eigen::Index dominance_active_count = 0;
        double additive_sum_squares = 0.0;
        double dominance_sum_squares = 0.0;
        for (const Eigen::Index marker : valid_indices)
        {
            const double old_additive = additive.coefficients(marker);
            const double old_dominance = dominance.coefficients(marker);
            const auto old_class_index
                = static_cast<std::size_t>(joint.assignment(marker));
            const double additive_linear
                = additive_projection.dot(marker, residual.adjusted_response)
                  + (additive_xtx_diag(marker) * old_additive);
            const double dominance_linear
                = dominance_projection.dot(marker, residual.adjusted_response)
                  + (dominance_xtx_diag(marker) * old_dominance);
            const auto sample = draw_marker(
                {.additive_quadratic = additive_xtx_diag(marker),
                 .additive_linear = additive_linear,
                 .dominance_quadratic = dominance_xtx_diag(marker),
                 .dominance_linear = dominance_linear,
                 .residual_variance = residual.variance},
                log_probabilities,
                joint.positive_probability,
                rng);
            const bool additive_active
                = sample.class_index == 1 || sample.class_index == 3;
            const bool dominance_active
                = sample.class_index == 2 || sample.class_index == 3;

            joint.dominance_sign(marker) = static_cast<std::uint8_t>(
                dominance_active && sample.positive_dominance);
            if constexpr (PositiveProbabilityUpdate == UpdatePolicy::Sampled)
            {
                if (dominance_active)
                {
                    positive_probability_updater_.observe(
                        sample.positive_dominance);
                }
            }

            additive.coefficients(marker) = sample.additive;
            dominance.coefficients(marker) = sample.dominance;
            joint.assignment(marker)
                = static_cast<std::uint8_t>(sample.class_index);
            apply_fitted_value_transition<GeneticMode::A, Layout>(
                additive_projection,
                marker,
                FittedValueTransition{
                    .old_class_index = old_class_index,
                    .new_class_index = sample.class_index,
                    .old_coefficient = old_additive,
                    .new_coefficient = sample.additive},
                additive.component_fitted_values,
                residual.adjusted_response);
            apply_fitted_value_transition<GeneticMode::D, Layout>(
                dominance_projection,
                marker,
                FittedValueTransition{
                    .old_class_index = old_class_index,
                    .new_class_index = sample.class_index,
                    .old_coefficient = old_dominance,
                    .new_coefficient = sample.dominance},
                dominance.component_fitted_values,
                residual.adjusted_response);

            if (additive_active)
            {
                ++additive_active_count;
                additive_sum_squares += sample.additive * sample.additive;
            }
            if (dominance_active)
            {
                ++dominance_active_count;
                dominance_sum_squares += sample.dominance * sample.dominance;
            }
            if constexpr (ProbabilitiesUpdate == UpdatePolicy::Sampled)
            {
                probabilities_updater_.observe(sample.class_index);
            }
        }

        additive_variance = additive_variance_sampler_(
            {.n = additive_active_count, .sum_squares = additive_sum_squares},
            rng);
        dominance_variance = dominance_variance_sampler_(
            {.n = dominance_active_count, .sum_squares = dominance_sum_squares},
            rng);
        if constexpr (ProbabilitiesUpdate == UpdatePolicy::Sampled)
        {
            probabilities_updater_.update(joint.probabilities, rng);
        }
        if constexpr (PositiveProbabilityUpdate == UpdatePolicy::Sampled)
        {
            positive_probability_updater_.update(
                joint.positive_probability, rng);
        }
    }

   private:
    auto draw_marker(
        const MarkerLikelihood& likelihood,
        const std::array<double, class_count>& log_probabilities,
        double positive_probability,
        std::mt19937_64& rng) -> MarkerSample
    {
        const auto additive_posterior = normal_.posterior_with_logL(
            {.quadratic = likelihood.additive_quadratic,
             .linear = likelihood.additive_linear,
             .scale = likelihood.residual_variance});
        const auto dominance_positive_posterior
            = half_normal_.posterior_with_logL(
                {.quadratic = likelihood.dominance_quadratic,
                 .linear = likelihood.dominance_linear,
                 .scale = likelihood.residual_variance},
                1);
        const auto dominance_negative_posterior
            = half_normal_.posterior_with_logL(
                {.quadratic = likelihood.dominance_quadratic,
                 .linear = likelihood.dominance_linear,
                 .scale = likelihood.residual_variance},
                -1);

        const double positive_log_weight
            = std::log(positive_probability)
              + dominance_positive_posterior.log_marginal_kernel;
        const double negative_log_weight
            = std::log1p(-positive_probability)
              + dominance_negative_posterior.log_marginal_kernel;
        const double dominance_log_likelihood
            = log_sum_exp(positive_log_weight, negative_log_weight);

        // NOIA A and D columns are orthogonal, so the AD marginal likelihood
        // factorizes into the two mode-local marginals.
        const std::array log_weights{
            log_probabilities[0],
            additive_posterior.log_likelihood_kernel + log_probabilities[1],
            dominance_log_likelihood + log_probabilities[2],
            additive_posterior.log_likelihood_kernel + dominance_log_likelihood
                + log_probabilities[3]};
        const std::size_t class_index = draw_class(log_weights, rng);
        const bool additive_active = class_index == 1 || class_index == 3;
        const bool dominance_active = class_index == 2 || class_index == 3;
        const bool positive_dominance
            = dominance_active
              && draw_positive_sign(
                  positive_log_weight, negative_log_weight, rng);
        const DominancePosterior& dominance_posterior
            = positive_dominance ? dominance_positive_posterior
                                 : dominance_negative_posterior;

        return {
            .class_index = class_index,
            .additive = additive_active
                            ? normal_.draw(additive_posterior.params, rng)
                            : 0.0,
            .dominance = dominance_active
                             ? half_normal_.draw(dominance_posterior, rng)
                             : 0.0,
            .positive_dominance = positive_dominance};
    }

    auto reset_samplers() -> void
    {
        additive_variance_sampler_.reset();
        dominance_variance_sampler_.reset();
        normal_.reset();
        half_normal_.reset();
        uniform_.reset();
    }

    [[nodiscard]] static auto log_sum_exp(double lhs, double rhs) noexcept
        -> double
    {
        const double maximum = std::max(lhs, rhs);
        return maximum
               + std::log(std::exp(lhs - maximum) + std::exp(rhs - maximum));
    }

    auto draw_class(
        std::array<double, class_count> log_weights,
        std::mt19937_64& rng) -> std::size_t
    {
        const double maximum = *std::ranges::max_element(log_weights);
        double total_weight = 0.0;
        for (double& weight : log_weights)
        {
            weight = std::exp(weight - maximum);
            total_weight += weight;
        }

        const double threshold = uniform_(rng) * total_weight;
        double cumulative_weight = 0.0;
        for (std::size_t class_index = 0; class_index < class_count;
             ++class_index)
        {
            cumulative_weight += log_weights[class_index];
            if (threshold < cumulative_weight)
            {
                return class_index;
            }
        }
        return class_count - 1;
    }

    auto draw_positive_sign(
        double positive_log_weight,
        double negative_log_weight,
        std::mt19937_64& rng) -> bool
    {
        const double maximum
            = std::max(positive_log_weight, negative_log_weight);
        const double positive_weight = std::exp(positive_log_weight - maximum);
        const double negative_weight = std::exp(negative_log_weight - maximum);
        return uniform_(rng) * (positive_weight + negative_weight)
               < positive_weight;
    }

    ScaledInvChi2Sampler<double> additive_variance_sampler_;
    ScaledInvChi2Sampler<double> dominance_variance_sampler_;
    [[no_unique_address]]
    SimplexUpdater<ProbabilitiesUpdate, class_count> probabilities_updater_;
    [[no_unique_address]]
    ProbabilityUpdater<PositiveProbabilityUpdate> positive_probability_updater_;
    NormalSampler<double> normal_{0.0};
    HalfNormalSampler<double> half_normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
};

template <
    UpdatePolicy ProbabilitiesUpdate,
    UpdatePolicy PositiveProbabilityUpdate>
[[nodiscard]] auto make_genetic_kernel(
    const JointTopology<
        GaussianPrior<VarianceLayout::Pooled>,
        JointSpikeSlabPrior<ProbabilitiesUpdate, PositiveProbabilityUpdate>>&
        prior)
    -> JointSpikeSlabKernel<ProbabilitiesUpdate, PositiveProbabilityUpdate>
{
    return JointSpikeSlabKernel<ProbabilitiesUpdate, PositiveProbabilityUpdate>{
        prior};
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_KERNEL_H_
