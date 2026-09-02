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
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <type_traits>

#include "gelex/bayes/detail/state_factory.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/probability_updater.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/half_normal_sampler.h"
#include "gelex/bayes/stats/normal_sampler.h"
#include "gelex/bayes/stats/scaled_inv_chi2_sampler.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

template <
    std::size_t ClassCount,
    UpdatePolicy ProbabilitiesUpdate,
    HalfNormalAsymmetry Axis>
class JointSpikeSlabKernel
{
    static_assert(ClassCount == 4);

    using AdditivePrior = GaussianPrior<VarianceLayout::Pooled>;
    using DominancePrior = HalfNormalPrior<Axis>;
    using ModePriors = ModeValues<
        GeneticMode::A | GeneticMode::D,
        AdditivePrior,
        DominancePrior>;
    using JointPrior = JointSpikeSlabPrior<ClassCount, ProbabilitiesUpdate>;
    using GeneticPrior = JointModeValues<ModePriors, JointPrior>;
    using GeneticState = genetic_state_t<GeneticPrior>;
    using JointState = JointSpikeSlabState<ClassCount>;
    using DominancePosterior = HalfNormalSampler<double>::Posterior;

    struct NoPositiveProbabilityUpdater
    {
    };

    using PositiveProbabilityUpdater = std::conditional_t<
        Axis == HalfNormalAsymmetry::Count,
        ProbabilityUpdater<UpdatePolicy::Sampled>,
        NoPositiveProbabilityUpdater>;

    static constexpr std::size_t class_count = ClassCount;
    static constexpr std::size_t negative_sign_index = 0;
    static constexpr std::size_t positive_sign_index = 1;

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

    struct CoefficientTransition
    {
        std::size_t old_class_index{};
        std::size_t new_class_index{};
        double old_coefficient{};
        double new_coefficient{};
    };

    struct DominanceParameters
    {
        std::array<double, 2> variances{};
        std::array<double, 2> sign_log_probabilities{};
    };

    struct EffectStatistics
    {
        Eigen::Index count{};
        double sum_squares{};

        auto observe(double coefficient) -> void
        {
            ++count;
            sum_squares += coefficient * coefficient;
        }
    };

    struct VarianceStatistics
    {
        EffectStatistics additive;
        EffectStatistics dominance;
        std::array<EffectStatistics, 2> dominance_by_sign{};
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
          positive_probability_updater_{
              make_positive_probability_updater(prior)},
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
        auto& additive_family = additive.family_state;
        auto& dominance_family = dominance.family_state;
        const auto& additive_xtx_diag = additive_projection.xtx_diag();
        const auto& dominance_xtx_diag = dominance_projection.xtx_diag();

        reset_updaters();

        std::array<double, class_count> log_probabilities{};
        for (std::size_t class_index = 0; class_index < class_count;
             ++class_index)
        {
            log_probabilities[class_index]
                = std::log(joint.probabilities[class_index]);
        }
        normal_.set_prior_var(additive_family.variance);
        const auto dominance_parameters
            = make_dominance_parameters(dominance_family);

        VarianceStatistics variance_statistics;
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
                dominance_parameters,
                rng);

            dominance_family.assignment(marker) = static_cast<std::uint8_t>(
                dominance_is_active(sample.class_index)
                    ? static_cast<std::uint8_t>(sample.positive_dominance) + 1
                    : 0);

            additive.coefficients(marker) = sample.additive;
            dominance.coefficients(marker) = sample.dominance;
            joint.assignment(marker)
                = static_cast<std::uint8_t>(sample.class_index);
            apply_fitted_value_transition<GeneticMode::A>(
                additive_projection,
                marker,
                {.old_class_index = old_class_index,
                 .new_class_index = sample.class_index,
                 .old_coefficient = old_additive,
                 .new_coefficient = sample.additive},
                additive_family.fitted_values,
                joint,
                residual.adjusted_response);
            apply_fitted_value_transition<GeneticMode::D>(
                dominance_projection,
                marker,
                {.old_class_index = old_class_index,
                 .new_class_index = sample.class_index,
                 .old_coefficient = old_dominance,
                 .new_coefficient = sample.dominance},
                dominance_family.fitted_values,
                joint,
                residual.adjusted_response);
            observe_sample(sample, variance_statistics);
        }

        update_variances(
            additive_family, dominance_family, variance_statistics, rng);
        update_probabilities(joint, dominance_family, rng);
    }

   private:
    [[nodiscard]] static auto make_positive_probability_updater(
        const GeneticPrior& prior) -> PositiveProbabilityUpdater
    {
        if constexpr (Axis == HalfNormalAsymmetry::Count)
        {
            return PositiveProbabilityUpdater{
                prior.mode_values()
                    .template get<GeneticMode::D>()
                    .positive_probability};
        }
        else
        {
            return {};
        }
    }

    [[nodiscard]] static auto make_dominance_parameters(
        const HalfNormalState<Axis>& state) -> DominanceParameters
    {
        if constexpr (Axis == HalfNormalAsymmetry::Count)
        {
            return {
                .variances = {state.variance, state.variance},
                .sign_log_probabilities
                = {std::log1p(-state.positive_probability),
                   std::log(state.positive_probability)}};
        }
        else
        {
            const double log_half = std::log(0.5);
            return {
                .variances = state.variances,
                .sign_log_probabilities = {log_half, log_half}};
        }
    }

    [[nodiscard]] static constexpr auto additive_is_active(
        std::size_t class_index) noexcept -> bool
    {
        return class_index == 1 || class_index == 3;
    }

    [[nodiscard]] static constexpr auto dominance_is_active(
        std::size_t class_index) noexcept -> bool
    {
        return class_index == 2 || class_index == 3;
    }

    auto draw_marker(
        const MarkerLikelihood& likelihood,
        const std::array<double, class_count>& log_probabilities,
        const DominanceParameters& dominance_parameters,
        std::mt19937_64& rng) -> MarkerSample
    {
        const auto additive_posterior = normal_.posterior_with_logL(
            {.quadratic = likelihood.additive_quadratic,
             .linear = likelihood.additive_linear,
             .scale = likelihood.residual_variance});
        const auto dominance_positive_posterior
            = half_normal_
                  .set_prior_var(
                      dominance_parameters.variances[positive_sign_index])
                  .posterior_with_logL(
                      {.quadratic = likelihood.dominance_quadratic,
                       .linear = likelihood.dominance_linear,
                       .scale = likelihood.residual_variance},
                      1);
        const auto dominance_negative_posterior
            = half_normal_
                  .set_prior_var(
                      dominance_parameters.variances[negative_sign_index])
                  .posterior_with_logL(
                      {.quadratic = likelihood.dominance_quadratic,
                       .linear = likelihood.dominance_linear,
                       .scale = likelihood.residual_variance},
                      -1);

        const double positive_log_weight
            = dominance_parameters.sign_log_probabilities[positive_sign_index]
              + dominance_positive_posterior.log_marginal_kernel;
        const double negative_log_weight
            = dominance_parameters.sign_log_probabilities[negative_sign_index]
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
        const bool additive_active = additive_is_active(class_index);
        const bool dominance_active = dominance_is_active(class_index);
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

    template <GeneticMode Mode>
        requires(Mode == GeneticMode::A || Mode == GeneticMode::D)
    [[nodiscard]] static constexpr auto fitted_component_index(
        std::size_t class_index) noexcept -> int
    {
        if (class_index >= class_count)
        {
            return JointState::no_component;
        }
        if constexpr (Mode == GeneticMode::A)
        {
            return JointState::additive_components[class_index];
        }
        else
        {
            return JointState::dominance_components[class_index];
        }
    }

    template <GeneticMode Mode>
        requires(Mode == GeneticMode::A || Mode == GeneticMode::D)
    static auto apply_fitted_value_transition(
        const bayes::GeneticProjection& projection,
        Eigen::Index marker,
        CoefficientTransition transition,
        Eigen::VectorXd& total_fitted_values,
        JointState& joint_state,
        Eigen::VectorXd& adjusted_response) -> void
    {
        std::array<bayes::AxpyTarget, 4> targets{};
        std::size_t target_count = 0;
        const auto append_target
            = [&](double scale, const Eigen::Ref<Eigen::VectorXd>& values)
        {
            assert(target_count < targets.size());
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            targets[target_count++] = bayes::AxpyTarget{scale, values};
        };
        const double coefficient_delta
            = transition.new_coefficient - transition.old_coefficient;
        if (coefficient_delta != 0.0)
        {
            append_target(-coefficient_delta, adjusted_response);
            append_target(coefficient_delta, total_fitted_values);
        }

        const auto append_joint_target
            = [&](std::size_t class_index, double delta)
        {
            const int component_index
                = fitted_component_index<Mode>(class_index);
            if (component_index == JointState::no_component || delta == 0.0)
            {
                return;
            }
            append_target(
                delta,
                joint_state.fitted_values.col(
                    static_cast<Eigen::Index>(component_index)));
        };

        const int old_component
            = fitted_component_index<Mode>(transition.old_class_index);
        const int new_component
            = fitted_component_index<Mode>(transition.new_class_index);
        if (old_component == new_component)
        {
            append_joint_target(transition.old_class_index, coefficient_delta);
        }
        else
        {
            append_joint_target(
                transition.old_class_index, -transition.old_coefficient);
            append_joint_target(
                transition.new_class_index, transition.new_coefficient);
        }

        if (target_count != 0)
        {
            projection.axpy(
                marker,
                std::span<const bayes::AxpyTarget>{
                    targets.data(), target_count});
        }
    }

    auto observe_sample(
        const MarkerSample& sample,
        VarianceStatistics& statistics) -> void
    {
        if (additive_is_active(sample.class_index))
        {
            statistics.additive.observe(sample.additive);
        }
        if (dominance_is_active(sample.class_index))
        {
            if constexpr (Axis == HalfNormalAsymmetry::Count)
            {
                statistics.dominance.observe(sample.dominance);
                positive_probability_updater_.observe(
                    sample.positive_dominance);
            }
            else
            {
                const std::size_t sign_index = sample.positive_dominance
                                                   ? positive_sign_index
                                                   : negative_sign_index;
                statistics.dominance_by_sign[sign_index].observe(
                    sample.dominance);
            }
        }
        if constexpr (ProbabilitiesUpdate == UpdatePolicy::Sampled)
        {
            probabilities_updater_.observe(sample.class_index);
        }
    }

    auto update_variances(
        GaussianState<VarianceLayout::Pooled>& additive_state,
        HalfNormalState<Axis>& dominance_state,
        const VarianceStatistics& statistics,
        std::mt19937_64& rng) -> void
    {
        additive_state.variance = additive_variance_sampler_(
            {.n = statistics.additive.count,
             .sum_squares = statistics.additive.sum_squares},
            rng);
        if constexpr (Axis == HalfNormalAsymmetry::Count)
        {
            dominance_state.variance = dominance_variance_sampler_(
                {.n = statistics.dominance.count,
                 .sum_squares = statistics.dominance.sum_squares},
                rng);
        }
        else
        {
            constexpr std::array sign_indices{
                negative_sign_index, positive_sign_index};
            for (const std::size_t sign_index : sign_indices)
            {
                const auto& sign_statistics
                    = statistics.dominance_by_sign[sign_index];
                dominance_state.variances[sign_index]
                    = dominance_variance_sampler_(
                        {.n = sign_statistics.count,
                         .sum_squares = sign_statistics.sum_squares},
                        rng);
            }
        }
    }

    auto update_probabilities(
        JointState& joint_state,
        HalfNormalState<Axis>& dominance_state,
        std::mt19937_64& rng) -> void
    {
        if constexpr (ProbabilitiesUpdate == UpdatePolicy::Sampled)
        {
            probabilities_updater_.update(joint_state.probabilities, rng);
        }
        if constexpr (Axis == HalfNormalAsymmetry::Count)
        {
            positive_probability_updater_.update(
                dominance_state.positive_probability, rng);
        }
    }

    auto reset_updaters() -> void
    {
        additive_variance_sampler_.reset();
        dominance_variance_sampler_.reset();
        normal_.reset();
        half_normal_.reset();
        uniform_.reset();
        if constexpr (ProbabilitiesUpdate == UpdatePolicy::Sampled)
        {
            probabilities_updater_.reset();
        }
        if constexpr (Axis == HalfNormalAsymmetry::Count)
        {
            positive_probability_updater_.reset();
        }
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
    PositiveProbabilityUpdater positive_probability_updater_;
    NormalSampler<double> normal_{0.0};
    HalfNormalSampler<double> half_normal_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_KERNEL_H_
