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

#include "gelex/model/bayes/prior.h"
#include <utility>

#include "gelex/exception.h"
#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/prior_config.h"
#include "gelex/model/bayes/prior_constants.h"

namespace gelex::bayes
{

namespace
{

auto compute_init_marker_variance(
    double target_variance,
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    double non_zero_proportion) -> double
{
    if (target_variance <= 0.0)
    {
        throw GelexException("Target variance must be positive");
    }
    if (non_zero_proportion <= 0.0)
    {
        throw GelexException("Non-zero marker proportion must be positive");
    }

    double total_genetic_variance = stats::detail::var(X).sum();
    auto num_non_zero_snps = total_genetic_variance * non_zero_proportion;

    if (num_non_zero_snps <= 0.0)
    {
        throw GelexException("Number of non-zero SNPs must be positive");
    }

    return target_variance / num_non_zero_snps;
}

auto make_variance_prior(double init) -> stats::detail::ScaledInvChiSqParams
{
    return {
        prior_constants::MARKER_VARIANCE_SHAPE,
        prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER * init};
}

}  // namespace

auto Priors::build_genetic_prior(
    const GeneticEffect& effect,
    const GeneticPriorConfig& prior,
    const PriorSetConfig& config) -> GeneticPrior
{
    const double target_variance
        = prior.heritability * config.phenotype_variance_;
    const auto& X = get_matrix_ref(effect.X);
    const auto num_cols = get_cols(effect.X);

    auto marker = [&]() -> MarkerPrior
    {
        switch (config.method_.base)
        {
            case BayesBase::A:
            case BayesBase::RR:
            {
                const double init = compute_init_marker_variance(
                    target_variance,
                    X,
                    prior_constants::NON_MIXTURE_PROPORTION);
                return ContinuousPrior{
                    {make_variance_prior(init), init, num_cols}};
            }
            case BayesBase::B:
            {
                const double non_zero = 1.0 - prior.proportion[0];
                const double init = compute_init_marker_variance(
                    target_variance, X, non_zero);
                return SpikePrior{
                    {make_variance_prior(init), init, num_cols},
                    {prior.proportion, config.method_.estimate_pi}};
            }
            case BayesBase::C:
            {
                const double non_zero = 1.0 - prior.proportion[0];
                const double init = compute_init_marker_variance(
                    target_variance, X, non_zero);
                return SpikePrior{
                    {make_variance_prior(init), init, 1},
                    {prior.proportion, config.method_.estimate_pi}};
            }
            case BayesBase::R:
            {
                const double non_zero = 1.0 - prior.proportion[0];
                const double init = compute_init_marker_variance(
                    target_variance, X, non_zero);
                return MixturePrior{
                    {make_variance_prior(init), init, 1},
                    {prior.proportion, true},
                    prior.multiplier};
            }
            case BayesBase::kCount:
                break;
        }
        std::unreachable();
    }();

    std::optional<SignPrior> sign;
    if (config.method_.dominance == bayes::DominancePolicy::asymmetric
        && effect.type == GeneticMode::D)
    {
        sign = SignPrior{config.positive_prob_};
    }

    return {effect.type, std::move(marker), sign};
}

Priors::Priors(
    const PriorSetConfig& config,
    const std::vector<GeneticEffect>& genetics_in,
    std::size_t num_random_effects)
    : residual_{
          {prior_constants::RESIDUAL_VARIANCE_SHAPE,
           prior_constants::RESIDUAL_VARIANCE_SCALE},
          config.residual_variance_proportion_ * config.phenotype_variance_,
          1}
{
    for (const auto& effect : genetics_in)
    {
        const auto* prior = config.find_genetic(effect.type);
        if (prior != nullptr)
        {
            genetics_.push_back(build_genetic_prior(effect, *prior, config));
        }
    }

    if (num_random_effects > 0)
    {
        const double target_variance
            = config.random_variance_proportion_ * config.phenotype_variance_;
        const double per_effect
            = target_variance / static_cast<double>(num_random_effects);

        for (std::size_t i = 0; i < num_random_effects; ++i)
        {
            random_.push_back(
                {{prior_constants::RANDOM_EFFECTS_SHAPE,
                  prior_constants::RANDOM_EFFECTS_SCALE},
                 per_effect,
                 1});
        }
    }
}

Priors::Priors(
    std::vector<GeneticPrior> genetics_in,
    std::vector<RandomPrior> random_in,
    ResidualPrior residual_in)
    : genetics_(std::move(genetics_in)),
      random_(std::move(random_in)),
      residual_(residual_in)
{
}

}  // namespace gelex::bayes
