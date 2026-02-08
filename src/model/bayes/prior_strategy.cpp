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

#include "gelex/model/bayes/prior_strategy.h"

#include "gelex/exception.h"
#include "gelex/model/bayes/prior_constants.h"
#include "utils/math_utils.h"

namespace gelex
{

PriorSetter::PriorSetter(PriorSpec spec) : spec_(spec) {}

auto PriorSetter::operator()(BayesModel& model, const PriorConfig& config)
    -> void
{
    if (model.additive() != nullptr)
    {
        apply_effect_prior(
            *model.additive(), spec_.additive, config.additive, config);
    }

    if (spec_.dominant.has_value() && model.dominant() != nullptr)
    {
        apply_effect_prior(
            *model.dominant(), *spec_.dominant, config.dominant, config);
    }

    const auto num_random_effect = static_cast<double>(model.random().size());
    if (!model.random().empty())
    {
        const double target_variance
            = config.random_variance_proportion * config.phenotype_variance;
        const double per_effect_variance = target_variance / num_random_effect;

        for (auto& effect : model.random())
        {
            set_random_effect_prior(per_effect_variance, effect);
        }
    }

    {
        const double target_variance
            = config.residual_variance_proportion * config.phenotype_variance;
        model.residual().prior
            = {prior_constants::RESIDUAL_VARIANCE_SHAPE,
               prior_constants::RESIDUAL_VARIANCE_SCALE};
        model.residual().init_variance = target_variance;
    }
}

auto PriorSetter::set_random_effect_prior(
    double variance,
    bayes::RandomEffect& effect) -> void
{
    effect.prior
        = {prior_constants::RANDOM_EFFECTS_SHAPE,
           prior_constants::RANDOM_EFFECTS_SCALE};
    effect.init_variance = variance;
}

auto compute_init_marker_variance(
    double target_variance,
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    double non_zero_marker_proption) -> double
{
    if (target_variance <= 0.0)
    {
        throw ArgumentValidationException("Target variance must be positive");
    }
    if (non_zero_marker_proption <= 0.0)
    {
        throw ArgumentValidationException(
            "Non-zero marker proportion must be positive");
    }

    double total_genetic_variance = gelex::detail::var(X).sum();
    auto num_non_zero_snps = total_genetic_variance * non_zero_marker_proption;

    if (num_non_zero_snps <= 0.0)
    {
        throw ArgumentValidationException(
            "Number of non-zero SNPs must be positive");
    }

    return target_variance / num_non_zero_snps;
}
}  // namespace gelex
