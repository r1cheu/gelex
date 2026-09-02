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

#ifndef GELEX_BAYES_GENETIC_KERNEL_HALF_NORMAL_UPDATER_H_
#define GELEX_BAYES_GENETIC_KERNEL_HALF_NORMAL_UPDATER_H_

#include <array>
#include <cstddef>
#include <random>

#include "gelex/bayes/genetic/kernel/allocation_updater.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/stats/half_normal_sampler.h"
#include "gelex/bayes/stats/scaled_inv_chi2_sampler.h"

namespace gelex::detail
{

class HalfNormalUpdater
{
    static constexpr std::size_t negative_index = 0;
    static constexpr std::size_t positive_index = 1;

    using CoefficientSampler = HalfNormalSampler<double>;
    using CoefficientPosterior = CoefficientSampler::Posterior;
    using SignAllocation
        = BinaryAllocationUpdater<MixtureWeightUpdate::Enabled>;
    using VarianceLikelihood = ScaledInvChi2Sampler<double>::Likelihood;

   public:
    struct Posterior
    {
        std::array<CoefficientPosterior, 2> coefficients;
        typename SignAllocation::Posterior sign;
        double log_marginal_kernel{};
    };

    struct Sample
    {
        double coefficient{};
        bool positive{};
    };

    explicit HalfNormalUpdater(const HalfNormalPrior& prior)
        : variance_sampler_{prior.variance.prior()},
          sign_allocation_{prior.positive_probability},
          coefficient_sampler_{prior.variance.initial_value()}
    {
    }

    auto begin_sweep(const HalfNormalState& state) -> void
    {
        variance_sampler_.reset();
        coefficient_sampler_.reset();
        statistics_ = {};
        variance_ = state.variance;
        sign_allocation_.begin_sweep(state.positive_probability);
    }

    [[nodiscard]] auto posterior_with_logL(
        const CoefficientSampler::Kernel& likelihood) -> Posterior
    {
        const auto negative_posterior
            = coefficient_sampler_.set_prior_var(variance_).posterior_with_logL(
                likelihood, -1);
        const auto positive_posterior
            = coefficient_sampler_.set_prior_var(variance_).posterior_with_logL(
                likelihood, 1);
        const auto sign_posterior = sign_allocation_.posterior(
            {negative_posterior.log_marginal_kernel,
             positive_posterior.log_marginal_kernel});

        return {
            .coefficients = {negative_posterior, positive_posterior},
            .sign = sign_posterior,
            .log_marginal_kernel = sign_posterior.log_marginal_kernel};
    }

    auto draw(const Posterior& posterior, std::mt19937_64& rng) -> Sample
    {
        const bool positive = sign_allocation_.draw(posterior.sign, rng);
        const std::size_t sign_index
            = positive ? positive_index : negative_index;
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
        const auto& coefficient_posterior = posterior.coefficients[sign_index];
        const double coefficient
            = coefficient_sampler_.draw(coefficient_posterior, rng);
        ++statistics_.n;
        statistics_.sum_squares += coefficient * coefficient;
        return {.coefficient = coefficient, .positive = positive};
    }

    auto update(HalfNormalState& state, std::mt19937_64& rng) -> void
    {
        state.variance = variance_sampler_(statistics_, rng);
        sign_allocation_.update(state.positive_probability, rng);
    }

   private:
    ScaledInvChi2Sampler<double> variance_sampler_;
    SignAllocation sign_allocation_;
    CoefficientSampler coefficient_sampler_;
    double variance_{};
    VarianceLikelihood statistics_{};
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_KERNEL_HALF_NORMAL_UPDATER_H_
