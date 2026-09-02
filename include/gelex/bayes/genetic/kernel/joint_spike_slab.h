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

#ifndef GELEX_BAYES_GENETIC_KERNEL_JOINT_SPIKE_SLAB_H_
#define GELEX_BAYES_GENETIC_KERNEL_JOINT_SPIKE_SLAB_H_

#include <Eigen/Core>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>

#include "gelex/bayes/detail/state_factory.h"
#include "gelex/bayes/genetic/kernel/allocation_updater.h"
#include "gelex/bayes/genetic/kernel/half_normal_updater.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/normal_sampler.h"
#include "gelex/bayes/stats/scaled_inv_chi2_sampler.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
class JointSpikeSlabKernel
{
    static_assert(ClassCount == 4);

    using AdditivePrior = GaussianPrior<VarianceLayout::Pooled>;
    using DominancePrior = HalfNormalPrior;
    using ModePriors = ModeValues<
        GeneticMode::A | GeneticMode::D,
        AdditivePrior,
        DominancePrior>;
    using JointPrior = JointSpikeSlabPrior<ClassCount, WeightUpdate>;
    using GeneticPrior = JointModeValues<ModePriors, JointPrior>;
    using GeneticState = genetic_state_t<GeneticPrior>;
    using JointState = JointSpikeSlabState<ClassCount>;
    using DominanceUpdater = HalfNormalUpdater;
    using DominanceSample = typename DominanceUpdater::Sample;

    static constexpr std::size_t class_count = ClassCount;

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
        DominanceSample dominance{};
    };

    struct CoefficientTransition
    {
        std::size_t old_class_index{};
        std::size_t new_class_index{};
        double old_coefficient{};
        double new_coefficient{};
    };

   public:
    explicit JointSpikeSlabKernel(const GeneticPrior& prior)
        : additive_variance_sampler_{prior.mode_values()
                                         .template get<GeneticMode::A>()
                                         .variance.prior()},
          allocation_updater_{prior.joint().probabilities},
          dominance_updater_{prior.mode_values().template get<GeneticMode::D>()}
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

        additive_variance_sampler_.reset();
        normal_.reset();
        allocation_updater_.begin_sweep(joint.probabilities);
        normal_.set_prior_var(additive_family.variance);
        dominance_updater_.begin_sweep(dominance_family);

        ScaledInvChi2Sampler<double>::Likelihood additive_statistics{};
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
                rng);

            dominance_family.assignment(marker) = static_cast<std::uint8_t>(
                is_active<GeneticMode::D>(sample.class_index)
                    ? static_cast<std::uint8_t>(sample.dominance.positive) + 1
                    : 0);

            additive.coefficients(marker) = sample.additive;
            dominance.coefficients(marker) = sample.dominance.coefficient;
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
                 .new_coefficient = sample.dominance.coefficient},
                dominance_family.fitted_values,
                joint,
                residual.adjusted_response);

            if (is_active<GeneticMode::A>(sample.class_index))
            {
                ++additive_statistics.n;
                additive_statistics.sum_squares
                    += sample.additive * sample.additive;
            }
        }

        additive_family.variance
            = additive_variance_sampler_(additive_statistics, rng);
        dominance_updater_.update(dominance_family, rng);
        allocation_updater_.update(joint.probabilities, rng);
    }

   private:
    template <GeneticMode Mode>
        requires(Mode == GeneticMode::A || Mode == GeneticMode::D)
    [[nodiscard]] static constexpr auto is_active(
        std::size_t class_index) noexcept -> bool
    {
        return fitted_component_index<Mode>(class_index)
               != JointState::no_component;
    }

    auto draw_marker(const MarkerLikelihood& likelihood, std::mt19937_64& rng)
        -> MarkerSample
    {
        const auto additive_posterior = normal_.posterior_with_logL(
            {.quadratic = likelihood.additive_quadratic,
             .linear = likelihood.additive_linear,
             .scale = likelihood.residual_variance});
        const auto dominance_posterior = dominance_updater_.posterior_with_logL(
            {.quadratic = likelihood.dominance_quadratic,
             .linear = likelihood.dominance_linear,
             .scale = likelihood.residual_variance});

        // NOIA A and D columns are orthogonal, so the AD marginal likelihood
        // factorizes into the two mode-local marginals.
        const auto allocation_posterior = allocation_updater_.posterior(
            std::array{
                0.0,
                additive_posterior.log_likelihood_kernel,
                dominance_posterior.log_marginal_kernel,
                additive_posterior.log_likelihood_kernel
                    + dominance_posterior.log_marginal_kernel});
        const std::size_t class_index
            = allocation_updater_.draw(allocation_posterior, rng);
        const bool additive_active = is_active<GeneticMode::A>(class_index);
        const bool dominance_active = is_active<GeneticMode::D>(class_index);
        const DominanceSample dominance_sample
            = dominance_active
                  ? dominance_updater_.draw(dominance_posterior, rng)
                  : DominanceSample{};

        return {
            .class_index = class_index,
            .additive = additive_active
                            ? normal_.draw(additive_posterior.params, rng)
                            : 0.0,
            .dominance = dominance_sample};
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

    ScaledInvChi2Sampler<double> additive_variance_sampler_;
    CategoricalAllocationUpdater<WeightUpdate, class_count> allocation_updater_;
    DominanceUpdater dominance_updater_;
    NormalSampler<double> normal_{0.0};
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_KERNEL_JOINT_SPIKE_SLAB_H_
