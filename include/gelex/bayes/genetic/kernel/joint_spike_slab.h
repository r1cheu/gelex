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
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <random>
#include <span>

#include "gelex/bayes/detail/normal_variance_conjugate_updater.h"
#include "gelex/bayes/detail/state_factory.h"
#include "gelex/bayes/genetic/kernel/coefficient_likelihood.h"
#include "gelex/bayes/genetic/kernel/dirichlet_conjugate_updater.h"
#include "gelex/bayes/genetic/kernel/probit_updater.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/half_quadratic_log_kernel.h"
#include "gelex/bayes/stats/log_categorical_distribution.h"
#include "gelex/bayes/stats/multi_quadratic_log_kernel.h"
#include "gelex/bayes/stats/quadratic_log_kernel.h"
#include "gelex/bayes/stats/truncated_normal_distribution.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/normal.h"

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
    using SignParameters = LogCategoricalDistribution<2>::param_type;

    static constexpr std::size_t class_count = ClassCount;
    static constexpr std::size_t negative_index = 0;
    static constexpr std::size_t positive_index = 1;

    struct DominancePosterior
    {
        HalfQuadraticLogKernel::Evaluation coefficient;
        SignParameters sign;
        double log_integral{};
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
        : additive_variance_updater_{prior.mode_values()
                                         .template get<GeneticMode::A>()
                                         .variance.prior},
          dominance_variance_updater_{prior.mode_values()
                                          .template get<GeneticMode::D>()
                                          .variance.prior},
          probit_updater_{make_multi_normal_prior(Eigen::Matrix2d::Identity())},
          probability_updater_{make_dirichlet_conjugate_updater<class_count>(
              prior.joint().probabilities)}
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
        assert(design.marker_covariate().has_value());
        const auto& marker_covariates = *design.marker_covariate();
        assert(marker_covariates.X().rows() == 2);
        auto& additive = state.mode_values().template get<GeneticMode::A>();
        auto& dominance = state.mode_values().template get<GeneticMode::D>();
        auto& joint = state.joint();
        auto& additive_family = additive.family_state;
        auto& dominance_family = dominance.family_state;
        std::normal_distribution<double> normal_distribution;
        TruncatedNormalDistribution<> dominance_distribution;
        TruncatedNormalDistribution<> probit_latent_distribution;

        const auto log_probabilities = make_log_weights(joint.probabilities);
        std::array<std::size_t, class_count> allocation_counts{};
        auto& probit_coefficients = dominance_family.probit_coefficients;

        const auto additive_prior = make_normal_prior(additive_family.variance);
        const auto dominance_prior
            = make_half_normal_prior(dominance_family.variance);

        std::size_t additive_count = 0;
        double additive_sum_squares = 0.0;
        std::size_t dominance_count = 0;
        double dominance_sum_squares = 0.0;
        Eigen::Matrix2d probit_likelihood_quadratic = Eigen::Matrix2d::Zero();
        Eigen::Vector2d probit_likelihood_linear = Eigen::Vector2d::Zero();
        for (const Eigen::Index marker : valid_indices)
        {
            const double old_additive = additive.coefficients(marker);
            const double old_dominance = dominance.coefficients(marker);
            const Eigen::Vector2d marker_covariate
                = marker_covariates.X().col(marker);
            const double linear_predictor
                = marker_covariate.dot(probit_coefficients);
            const std::array<double, 2> dominance_log_probabilities{
                log_norm_cdf(-linear_predictor),
                log_norm_cdf(linear_predictor)};
            const auto old_class_index
                = static_cast<std::size_t>(joint.assignment(marker));
            const auto additive_likelihood = make_coefficient_likelihood(
                additive_projection, marker, old_additive, residual);
            const auto dominance_likelihood = make_coefficient_likelihood(
                dominance_projection, marker, old_dominance, residual);
            const auto additive_posterior
                = additive_likelihood + additive_prior;
            const double additive_log_integral
                = additive_posterior.log_integral();
            const auto dominance_posterior = make_dominance_posterior(
                dominance_likelihood,
                dominance_prior,
                dominance_log_probabilities);

            // NOIA A and D columns are orthogonal, so the AD marginal
            // likelihood factorizes into the two mode-local marginals.
            const auto allocation_parameters = make_mixture_posterior_weights(
                log_probabilities,
                std::array{
                    0.0,
                    additive_log_integral,
                    dominance_posterior.log_integral,
                    additive_log_integral + dominance_posterior.log_integral});
            const std::size_t class_index
                = allocation_distribution_(rng, allocation_parameters);
            if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
            {
                // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
                ++allocation_counts[class_index];
            }

            const bool additive_active = is_active<GeneticMode::A>(class_index);
            const bool dominance_active
                = is_active<GeneticMode::D>(class_index);
            double new_additive = 0.0;
            double new_dominance = 0.0;
            if (additive_active)
            {
                new_additive = normal_distribution(
                    rng, additive_posterior.normal_parameters());
                ++additive_count;
                additive_sum_squares += new_additive * new_additive;
            }
            if (dominance_active)
            {
                const std::size_t sign_index
                    = sign_distribution_(rng, dominance_posterior.sign);
                const HalfLine support = sign_index == positive_index
                                             ? HalfLine::Positive
                                             : HalfLine::Negative;
                new_dominance = dominance_distribution(
                    rng,
                    dominance_posterior.coefficient.truncated_normal_parameters(
                        support));
                const double latent_value = probit_latent_distribution(
                    rng,
                    TruncatedNormalDistribution<>::param_type{
                        linear_predictor, 1.0, support});
                probit_likelihood_quadratic.noalias()
                    += marker_covariate * marker_covariate.transpose();
                probit_likelihood_linear.noalias()
                    += marker_covariate * latent_value;
                ++dominance_count;
                dominance_sum_squares += new_dominance * new_dominance;
            }

            additive.coefficients(marker) = new_additive;
            dominance.coefficients(marker) = new_dominance;
            joint.assignment(marker) = static_cast<std::uint8_t>(class_index);
            apply_fitted_value_transition<GeneticMode::A>(
                additive_projection,
                marker,
                {.old_class_index = old_class_index,
                 .new_class_index = class_index,
                 .old_coefficient = old_additive,
                 .new_coefficient = new_additive},
                additive_family.fitted_values,
                joint,
                residual.adjusted_response);
            apply_fitted_value_transition<GeneticMode::D>(
                dominance_projection,
                marker,
                {.old_class_index = old_class_index,
                 .new_class_index = class_index,
                 .old_coefficient = old_dominance,
                 .new_coefficient = new_dominance},
                dominance_family.fitted_values,
                joint,
                residual.adjusted_response);
        }

        if (dominance_count != 0)
        {
            probit_updater_.update(
                probit_coefficients,
                MultiQuadraticLogKernel{
                    probit_likelihood_quadratic, probit_likelihood_linear, 0.0},
                rng);
        }
        additive_variance_updater_.update(
            additive_family.variance,
            additive_count,
            additive_sum_squares,
            rng);
        dominance_variance_updater_.update(
            dominance_family.variance,
            dominance_count,
            dominance_sum_squares,
            rng);

        probability_updater_.update(
            joint.probabilities, allocation_counts, rng);
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

    [[nodiscard]] static auto make_dominance_posterior(
        const QuadraticLogKernel& likelihood,
        const HalfQuadraticLogKernel& prior,
        const std::array<double, 2>& log_probabilities) -> DominancePosterior
    {
        const auto coefficient_posterior = (likelihood + prior).evaluate();
        const std::array component_log_integrals{
            coefficient_posterior.log_integral(HalfLine::Negative),
            coefficient_posterior.log_integral(HalfLine::Positive)};
        const auto sign_parameters = make_mixture_posterior_weights(
            log_probabilities, component_log_integrals);

        const double negative_log_weight
            = log_probabilities[negative_index]
              + component_log_integrals[negative_index];
        const double positive_log_weight
            = log_probabilities[positive_index]
              + component_log_integrals[positive_index];
        const double maximum_log_weight
            = std::max(negative_log_weight, positive_log_weight);
        const double log_integral
            = maximum_log_weight
              + std::log(
                  std::exp(negative_log_weight - maximum_log_weight)
                  + std::exp(positive_log_weight - maximum_log_weight));
        return {
            .coefficient = coefficient_posterior,
            .sign = sign_parameters,
            .log_integral = log_integral};
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

    NormalVarianceConjugateUpdater additive_variance_updater_;
    NormalVarianceConjugateUpdater dominance_variance_updater_;
    ProbitUpdater probit_updater_;
    [[no_unique_address]] DirichletConjugateUpdater<class_count, WeightUpdate>
        probability_updater_;
    LogCategoricalDistribution<class_count> allocation_distribution_;
    LogCategoricalDistribution<2> sign_distribution_;
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_KERNEL_JOINT_SPIKE_SLAB_H_
