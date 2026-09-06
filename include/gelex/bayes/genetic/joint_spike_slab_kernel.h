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

#include "gelex/bayes/basic_state.h"
#include "gelex/bayes/detail/normal_variance_conjugate_updater.h"
#include "gelex/bayes/detail/state_factory.h"
#include "gelex/bayes/genetic/detail/apply_fitted_update.h"
#include "gelex/bayes/genetic/detail/coefficient_likelihood.h"
#include "gelex/bayes/genetic/detail/dirichlet_conjugate_updater.h"
#include "gelex/bayes/genetic/detail/probit_updater.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic_policy.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/operations.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/stats/half_quadratic_log_kernel.h"
#include "gelex/bayes/stats/log_categorical_distribution.h"
#include "gelex/bayes/stats/multi_quadratic_log_kernel.h"
#include "gelex/bayes/stats/quadratic_log_kernel.h"
#include "gelex/bayes/stats/truncated_normal_distribution.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/normal.h"

namespace gelex::detail
{

template <MixtureWeightUpdate WeightUpdate>
class JointSpikeSlabKernel
{
    using AdditivePrior = GaussianPrior<VarianceLayout::Pooled>;
    using DominancePrior = HalfNormalPrior;
    using ModePriors = ModeValues<
        GeneticMode::A | GeneticMode::D,
        AdditivePrior,
        DominancePrior>;
    using JointPrior = JointSpikeSlabPrior<WeightUpdate>;
    using GeneticPrior = JointModeValues<ModePriors, JointPrior>;
    using GeneticState = genetic_state_t<GeneticPrior>;
    using JointState = JointSpikeSlabState;
    using SignParameters = LogCategoricalDistribution<2>::param_type;

    static constexpr std::size_t class_count
        = JointSpikeSlabSpec<>::class_count;
    static constexpr std::size_t negative_index = 0;
    static constexpr std::size_t positive_index = 1;

    struct DominancePosterior
    {
        HalfQuadraticLogKernel::Evaluation coefficient;
        SignParameters sign;
        double log_integral{};
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
        auto& additive = state.template get<GeneticMode::A>();
        auto& dominance = state.template get<GeneticMode::D>();
        auto& joint = state.joint();
        std::normal_distribution<double> normal_distribution;
        TruncatedNormalDistribution<> dominance_distribution;
        TruncatedNormalDistribution<> probit_latent_distribution;

        const auto log_probabilities = make_log_weights(joint.probabilities());
        auto& probit_coefficients = dominance.probit_coefficients();

        const auto additive_prior = make_normal_prior(additive.variance());
        const auto dominance_prior
            = make_half_normal_prior(dominance.variance());

        double additive_sum_squares = 0.0;
        double dominance_sum_squares = 0.0;
        Eigen::Matrix2d probit_likelihood_quadratic = Eigen::Matrix2d::Zero();
        Eigen::Vector2d probit_likelihood_linear = Eigen::Vector2d::Zero();
        additive_fitted_delta_.setZero(design.rows());
        dominance_fitted_delta_.setZero(design.rows());
        for (const Eigen::Index marker : valid_indices)
        {
            const double old_additive = additive.coefficients()(marker);
            const double old_dominance = dominance.coefficients()(marker);
            const Eigen::Vector2d marker_covariate
                = marker_covariates.X().col(marker);
            const double linear_predictor
                = marker_covariate.dot(probit_coefficients);
            const std::array<double, 2> dominance_log_probabilities{
                log_norm_cdf(-linear_predictor),
                log_norm_cdf(linear_predictor)};
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
            const bool additive_active = is_active<GeneticMode::A>(class_index);
            const bool dominance_active
                = is_active<GeneticMode::D>(class_index);
            double new_additive = 0.0;
            double new_dominance = 0.0;
            if (additive_active)
            {
                new_additive = normal_distribution(
                    rng, additive_posterior.normal_parameters());
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
                dominance_sum_squares += new_dominance * new_dominance;
            }

            const auto updates = joint.transition(
                marker,
                static_cast<std::uint8_t>(class_index),
                typename JointState::ModeCoefficients{
                    old_additive, old_dominance},
                typename JointState::ModeCoefficients{
                    new_additive, new_dominance});
            additive.transition(marker, new_additive);
            dominance.transition(marker, new_dominance);
            const double additive_delta = new_additive - old_additive;
            const std::array additive_targets{
                bayes::AxpyTarget{-additive_delta, residual.adjusted_response},
                bayes::AxpyTarget{additive_delta, additive_fitted_delta_}};
            apply_fitted_update(
                additive_projection,
                marker,
                updates.template get<GeneticMode::A>(),
                std::span{additive_targets});
            const double dominance_delta = new_dominance - old_dominance;
            const std::array dominance_targets{
                bayes::AxpyTarget{-dominance_delta, residual.adjusted_response},
                bayes::AxpyTarget{dominance_delta, dominance_fitted_delta_}};
            apply_fitted_update(
                dominance_projection,
                marker,
                updates.template get<GeneticMode::D>(),
                std::span{dominance_targets});
        }

        additive.transition(additive_fitted_delta_);
        dominance.transition(dominance_fitted_delta_);
        const auto& counts = joint.class_counts();
        const auto additive_count = counts[1] + counts[3];
        const auto dominance_count = counts[2] + counts[3];
        if (dominance_count != 0)
        {
            probit_updater_.update(
                probit_coefficients,
                MultiQuadraticLogKernel{
                    probit_likelihood_quadratic, probit_likelihood_linear, 0.0},
                rng);
        }
        additive_variance_updater_.update(
            additive.variance(), additive_count, additive_sum_squares, rng);
        dominance_variance_updater_.update(
            dominance.variance(), dominance_count, dominance_sum_squares, rng);

        auto allocation_counts = counts;
        // Skipped markers remain NULL but do not contribute to the posterior.
        allocation_counts[0]
            -= static_cast<std::size_t>(joint.assignments().size())
               - valid_indices.size();
        probability_updater_.update(
            joint.probabilities(), allocation_counts, rng);
    }

   private:
    template <GeneticMode Mode>
        requires(Mode == GeneticMode::A || Mode == GeneticMode::D)
    [[nodiscard]] static constexpr auto is_active(
        std::size_t class_index) noexcept -> bool
    {
        return JointState::template fitted_component_index<Mode>(class_index)
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

    Eigen::VectorXd additive_fitted_delta_;
    Eigen::VectorXd dominance_fitted_delta_;
    NormalVarianceConjugateUpdater additive_variance_updater_;
    NormalVarianceConjugateUpdater dominance_variance_updater_;
    ProbitUpdater probit_updater_;
    [[no_unique_address]] DirichletConjugateUpdater<class_count, WeightUpdate>
        probability_updater_;
    LogCategoricalDistribution<class_count> allocation_distribution_;
    LogCategoricalDistribution<2> sign_distribution_;
};

template <MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_kernel(
    const JointModeValues<
        ModeValues<
            GeneticMode::A | GeneticMode::D,
            GaussianPrior<VarianceLayout::Pooled>,
            HalfNormalPrior>,
        JointSpikeSlabPrior<WeightUpdate>>& prior)
{
    return JointSpikeSlabKernel<WeightUpdate>{prior};
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_KERNEL_H_
