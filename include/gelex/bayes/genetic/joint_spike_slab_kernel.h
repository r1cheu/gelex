// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
#include <type_traits>
#include <variant>

#include "gelex/bayes/detail/normal_variance_conjugate_updater.h"
#include "gelex/bayes/genetic/detail/coefficient_likelihood.h"
#include "gelex/bayes/genetic/detail/dirichlet_conjugate_updater.h"
#include "gelex/bayes/genetic/detail/probit_updater.h"
#include "gelex/bayes/genetic/factory.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/half_quadratic_log_kernel.h"
#include "gelex/bayes/stats/log_categorical_distribution.h"
#include "gelex/bayes/stats/multi_quadratic_log_kernel.h"
#include "gelex/bayes/stats/quadratic_log_kernel.h"
#include "gelex/bayes/stats/truncated_normal_distribution.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/normal.h"

namespace gelex
{

template <MixtureWeightUpdate WeightUpdate>
class JointSpikeSlabKernel
{
    using additive_prior_type = GaussianPrior<VarianceLayout::Pooled>;
    using dominance_prior_type = HalfNormalPrior;
    using mode_priors_type = ModeValues<
        GeneticMode::A | GeneticMode::D,
        additive_prior_type,
        dominance_prior_type>;
    using joint_prior_type = JointSpikeSlabPrior<WeightUpdate>;
    using genetic_prior_type
        = JointModeValues<mode_priors_type, joint_prior_type>;
    using genetic_state_type = genetic_state_t<genetic_prior_type>;
    using sign_parameters_type = LogCategoricalDistribution<2>::param_type;

    static constexpr std::size_t class_count
        = JointSpikeSlabSpec<>::class_count;
    using probability_updater_type = std::conditional_t<
        WeightUpdate == MixtureWeightUpdate::Enabled,
        detail::DirichletConjugateUpdater<class_count>,
        std::monostate>;

    static constexpr std::size_t negative_index = 0;
    static constexpr std::size_t positive_index = 1;

    struct DominancePosterior
    {
        HalfQuadraticLogKernel::Evaluation coefficient;
        sign_parameters_type sign;
        double log_integral{};
    };

   public:
    explicit JointSpikeSlabKernel(const genetic_prior_type& prior)
        : additive_variance_updater_{prior.mode_values()
                                         .template get<GeneticMode::A>()
                                         .variance.prior},
          dominance_variance_updater_{prior.mode_values()
                                          .template get<GeneticMode::D>()
                                          .variance.prior},
          probit_updater_{make_multi_normal_prior(Eigen::Matrix2d::Identity())},
          probability_updater_{
              [&]() -> probability_updater_type
              {
                  if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
                  {
                      return probability_updater_type{
                          prior.joint().probabilities.prior};
                  }
                  else
                  {
                      return {};
                  }
              }()}
    {
    }

    auto step(
        const bayes::GeneticDesign& design,
        genetic_state_type& state,
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
        auto& annotation_coefficients = dominance.annotation_coefficients();

        const auto additive_prior = make_normal_prior(additive.variance());
        const auto dominance_prior
            = make_half_normal_prior(dominance.variance());

        double additive_sum_squares = 0.0;
        double dominance_sum_squares = 0.0;
        Eigen::Matrix2d probit_likelihood_quadratic = Eigen::Matrix2d::Zero();
        Eigen::Vector2d probit_likelihood_linear = Eigen::Vector2d::Zero();
        for (const Eigen::Index marker : valid_indices)
        {
            const double old_additive = additive.coefficients()(marker);
            const double old_dominance = dominance.coefficients()(marker);
            const Eigen::Vector2d marker_covariate
                = marker_covariates.X().col(marker);
            const double linear_predictor
                = marker_covariate.dot(annotation_coefficients);
            const std::array<double, 2> dominance_log_probabilities{
                log_norm_cdf(-linear_predictor),
                log_norm_cdf(linear_predictor)};
            const auto additive_likelihood
                = detail::make_coefficient_likelihood(
                    additive_projection, marker, old_additive, residual);
            const auto dominance_likelihood
                = detail::make_coefficient_likelihood(
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

            if (old_additive != new_additive)
            {
                additive_projection.axpy(
                    marker,
                    old_additive - new_additive,
                    residual.adjusted_response);
            }
            if (old_dominance != new_dominance)
            {
                dominance_projection.axpy(
                    marker,
                    old_dominance - new_dominance,
                    residual.adjusted_response);
            }
            additive.transition(marker, new_additive);
            dominance.transition(marker, new_dominance);
            joint.transition(marker, static_cast<std::uint8_t>(class_index));
        }

        const auto& counts = joint.class_counts();
        const auto additive_count = counts[1] + counts[3];
        const auto dominance_count = counts[2] + counts[3];
        if (dominance_count != 0)
        {
            probit_updater_.update(
                annotation_coefficients,
                MultiQuadraticLogKernel{
                    probit_likelihood_quadratic, probit_likelihood_linear, 0.0},
                rng);
        }
        additive_variance_updater_.update(
            additive.variance(), additive_count, additive_sum_squares, rng);
        dominance_variance_updater_.update(
            dominance.variance(), dominance_count, dominance_sum_squares, rng);

        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            auto allocation_counts = counts;
            // Skipped markers remain NULL but do not contribute to the
            // posterior.
            allocation_counts[0]
                -= static_cast<std::size_t>(joint.assignments().size())
                   - valid_indices.size();
            joint.set_probabilities(
                probability_updater_.draw(allocation_counts, rng));
        }
    }

   private:
    template <GeneticMode Mode>
        requires(Mode == GeneticMode::A || Mode == GeneticMode::D)
    [[nodiscard]] static constexpr auto is_active(
        std::size_t class_index) noexcept -> bool
    {
        return joint_class_activates(class_index, Mode);
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

    detail::NormalVarianceConjugateUpdater additive_variance_updater_;
    detail::NormalVarianceConjugateUpdater dominance_variance_updater_;
    detail::ProbitUpdater probit_updater_;
    [[no_unique_address]] probability_updater_type probability_updater_;
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

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_JOINT_SPIKE_SLAB_KERNEL_H_
