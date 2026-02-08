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
#include <optional>

#include "../src/types/bayes_effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior_constants.h"

namespace gelex
{

struct Prior
{
    Eigen::VectorXd mixture_proportions;
    Eigen::VectorXd mixture_scales;
    double heritability;
};

auto compute_init_marker_variance(
    double target_variance,
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    double non_zero_marker_proption) -> double;

struct PriorConfig
{
    double phenotype_variance{0.0};
    Prior additive{Eigen::VectorXd::Zero(2), Eigen::VectorXd(5), 0.5};
    Prior dominant{Eigen::VectorXd::Zero(2), Eigen::VectorXd(5), 0.2};
    double random_variance_proportion{0.1};
    double residual_variance_proportion{0.3};
};

enum class PriorType : uint8_t
{
    NonMixture,    // BayesA, BayesRR: all markers non-zero
    PiMixture,     // BayesB, BayesC: mixture with pi
    ScaleMixture,  // BayesR: mixture with scaled variances
};

enum class VarianceScope : uint8_t
{
    PerMarker,  // BayesA, BayesB: individual marker variances
    Shared,     // BayesRR, BayesC, BayesR: shared variance
};

struct EffectPriorSpec
{
    PriorType type;
    VarianceScope scope;
    bool estimate_pi;
};

struct PriorSpec
{
    EffectPriorSpec additive;
    std::optional<EffectPriorSpec> dominant;
};

class PriorSetter
{
   public:
    explicit PriorSetter(PriorSpec spec);
    auto operator()(BayesModel& model, const PriorConfig& config) -> void;

   private:
    PriorSpec spec_;

    template <typename EffectT>
    static auto apply_effect_prior(
        EffectT& effect,
        const EffectPriorSpec& spec,
        const Prior& effect_prior,
        const PriorConfig& config) -> void;

    static auto set_random_effect_prior(
        double variance,
        bayes::RandomEffect& effect) -> void;
};

// --- template implementation ---

template <typename EffectT>
auto PriorSetter::apply_effect_prior(
    EffectT& effect,
    const EffectPriorSpec& spec,
    const Prior& effect_prior,
    const PriorConfig& config) -> void
{
    const double target_variance
        = effect_prior.heritability * config.phenotype_variance;

    switch (spec.type)
    {
        case PriorType::NonMixture:
        {
            const double init_marker_variance = compute_init_marker_variance(
                target_variance,
                bayes::get_matrix_ref(effect.X),
                prior_constants::NON_MIXTURE_PROPORTION);

            effect.init_marker_variance = init_marker_variance;
            effect.marker_variance_prior
                = {prior_constants::MARKER_VARIANCE_SHAPE,
                   prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
                       * init_marker_variance};
            effect.marker_variance_size = bayes::get_cols(effect.X);
            break;
        }
        case PriorType::PiMixture:
        {
            const double non_mixture_prop
                = 1.0 - effect_prior.mixture_proportions[0];
            const double init_marker_variance = compute_init_marker_variance(
                target_variance,
                bayes::get_matrix_ref(effect.X),
                non_mixture_prop);

            effect.init_marker_variance = init_marker_variance;
            effect.marker_variance_prior
                = {prior_constants::MARKER_VARIANCE_SHAPE,
                   prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
                       * init_marker_variance};
            effect.init_pi.emplace(effect_prior.mixture_proportions);
            effect.marker_variance_size
                = (spec.scope == VarianceScope::PerMarker)
                      ? bayes::get_cols(effect.X)
                      : 1;
            break;
        }
        case PriorType::ScaleMixture:
        {
            const double non_mixture_prop
                = 1.0 - effect_prior.mixture_proportions[0];
            const double init_marker_variance = compute_init_marker_variance(
                target_variance,
                bayes::get_matrix_ref(effect.X),
                non_mixture_prop);

            effect.init_marker_variance = init_marker_variance;
            effect.marker_variance_prior
                = {prior_constants::MARKER_VARIANCE_SHAPE,
                   prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER
                       * init_marker_variance};
            effect.init_pi.emplace(effect_prior.mixture_proportions);
            effect.scale.emplace(effect_prior.mixture_scales);
            effect.marker_variance_size = 1;
            break;
        }
    }

    if (spec.estimate_pi)
    {
        effect.estimate_pi = true;
    }
}

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_PRIOR_STRATEGY_H_
