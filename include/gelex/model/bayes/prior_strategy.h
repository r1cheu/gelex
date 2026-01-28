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

#ifndef GELEX_MODEL_BAYES_PRIOR_STRATEGY_H_
#define GELEX_MODEL_BAYES_PRIOR_STRATEGY_H_

#include <Eigen/Core>

#include "../src/types/bayes_effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior_constants.h"

namespace gelex
{

/**
 * Configuration for a single effect type prior
 */
struct Prior
{
    Eigen::VectorXd
        mixture_proportions;  ///< Mixture proportions for effect sizes
    Eigen::VectorXd
        mixture_scales;   ///< Scale parameters for mixture components
    double heritability;  ///< Heritability proportion for this effect type
};

auto compute_init_marker_variance(
    double target_variance,
    const Eigen::Ref<const Eigen::MatrixXd>& design_matrix,
    double non_zero_marker_proption) -> double;
// use for unstandardized genotype matrix
auto compute_init_marker_variance(
    double target_variance,
    const Eigen::Ref<const Eigen::VectorXd>& genetic_variance,
    double non_zero_marker_proption) -> double;

/**
 * Complete prior configuration for Bayesian models
 *
 * Contains prior settings for all effect types and variance components
 */
struct PriorConfig
{
    double phenotype_variance{0.0};  ///< Total phenotypic variance
    Prior additive{
        Eigen::VectorXd::Zero(2),
        Eigen::VectorXd(5),
        0.5};  ///< Additive effect prior configuration
    Prior dominant{
        Eigen::VectorXd::Zero(2),
        Eigen::VectorXd(5),
        0.2};  ///< Dominant effect prior configuration
    double random_variance_proportion{
        0.1};  ///< Proportion of variance for random effects
    double residual_variance_proportion{
        0.3};  ///< Proportion of variance for residuals
};

/**
 * Abstract base class for setting priors in Bayesian models
 *
 * Implements the Strategy pattern for different prior configurations
 */
class PriorSetter
{
   public:
    PriorSetter() = default;
    PriorSetter(const PriorSetter&) = default;
    PriorSetter(PriorSetter&&) = default;
    PriorSetter& operator=(const PriorSetter&) = default;
    PriorSetter& operator=(PriorSetter&&) = default;
    virtual ~PriorSetter() = default;

    virtual auto operator()(BayesModel& model, const PriorConfig& config)
        -> void;

   protected:
    static auto set_random_effect_prior(
        double variance,
        bayes::RandomEffect& effect) -> void;

   private:
    virtual auto set_additive_effect_prior(
        bayes::AdditiveEffect& effect,
        const PriorConfig& prior) -> void
        = 0;
    virtual auto set_dominant_effect_prior(
        bayes::DominantEffect& /*effect*/,
        const PriorConfig& /*prior*/) -> void {
        // Default implementation does nothing
    };
};
}  // namespace gelex
namespace gelex::detail
{

template <typename T>
concept CppBayesEffectCore = requires(T effect, double var, size_t size) {
    { effect.design_matrix };
    { effect.init_marker_variance = var };
    { effect.marker_variance_size = size };
    { effect.marker_variance_prior };
};

template <typename T>
concept CppBayesEffectWithPi = CppBayesEffectCore<T> && requires(T effect) {
    { effect.init_pi };
};

template <typename T>
concept CppBayesEffectWithScale = CppBayesEffectCore<T> && requires(T effect) {
    { effect.scale };
};

template <CppBayesEffectCore EffectT>
auto set_core_effect_prior(
    EffectT& effect,
    double target_variance,
    double non_mixture_proportion) -> void
{
    const double init_marker_variance = gelex::compute_init_marker_variance(
        target_variance,
        gelex::bayes::get_matrix_ref(effect.design_matrix),
        non_mixture_proportion);

    effect.init_marker_variance = init_marker_variance;
    effect.marker_variance_prior
        = {gelex::prior_constants::MARKER_VARIANCE_SHAPE,
           gelex::prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
               * init_marker_variance};
}

template <CppBayesEffectCore EffectT, typename EffectPriorT>
auto set_non_mixture_prior_strategy(
    EffectT& effect,
    const EffectPriorT& effect_prior,
    const gelex::PriorConfig& prior) -> void
{
    const double target_variance
        = effect_prior.heritability * prior.phenotype_variance;

    set_core_effect_prior(
        effect,
        target_variance,
        gelex::prior_constants::NON_MIXTURE_PROPORTION);

    effect.marker_variance_size = gelex::bayes::get_cols(effect.design_matrix);
}

template <CppBayesEffectWithPi EffectT, typename EffectPriorT>
auto set_pi_mixture_prior_strategy(
    EffectT& effect,
    const EffectPriorT& effect_prior,
    const gelex::PriorConfig& prior,
    size_t marker_variance_size) -> void
{
    const double target_variance
        = effect_prior.heritability * prior.phenotype_variance;
    const double non_mixture_prop = 1.0 - effect_prior.mixture_proportions[0];

    set_core_effect_prior(effect, target_variance, non_mixture_prop);

    effect.init_pi.emplace(effect_prior.mixture_proportions);
    effect.marker_variance_size = marker_variance_size;
}

template <CppBayesEffectWithScale EffectT, typename EffectPriorT>
auto set_scale_mixture_prior_strategy(
    EffectT& effect,
    const EffectPriorT& effect_prior,
    const gelex::PriorConfig& prior) -> void
{
    const double target_variance
        = effect_prior.heritability * prior.phenotype_variance;
    const double non_mixture_prop = 1.0 - effect_prior.mixture_proportions[0];

    set_core_effect_prior(effect, target_variance, non_mixture_prop);

    effect.init_pi.emplace(effect_prior.mixture_proportions);
    effect.scale.emplace(effect_prior.mixture_scales);
    effect.marker_variance_size = 1;
}

}  // namespace gelex::detail

#endif  // GELEX_MODEL_BAYES_PRIOR_STRATEGY_H_
