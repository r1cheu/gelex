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

#ifndef GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_VARIANCE_H_
#define GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_VARIANCE_H_

#include <array>
#include <random>

#include "gelex/algo/infer/mcmc/kernels/common.h"
#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/exception.h"
#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class SingleGeneticVarianceStep final : public Step
{
   public:
    SingleGeneticVarianceStep(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        std::mt19937_64& rng)
        : block_(block),
          marker_variance_(
              prior.query<bayes::SingleMarkerVarianceCap>() != nullptr
                  ? prior.query<bayes::SingleMarkerVarianceCap>()->variance()
                  : throw GelexException(
                        "SingleGeneticVarianceStep: prior lacks marker "
                        "variance capability")),
          variance_(block.prior_state()
                        .require<bayes::SingleVarianceStateCap>()
                        .variance()),
          proportion_cap_(
              block.prior_state().query<bayes::SingleProportionStateCap>()),
          multiplier_cap_(prior.query<bayes::SingleMultiplierCap>()),
          rng_(rng)
    {
    }

    auto step() -> void override
    {
        auto sampler = make_variance_sampler(marker_variance_.parameter());
        const auto& coeffs = block_.state().coeffs;
        auto* proportion = proportion_cap_ == nullptr
                               ? nullptr
                               : &proportion_cap_->proportion();
        const auto* multiplier = multiplier_cap_ == nullptr
                                     ? nullptr
                                     : &multiplier_cap_->multiplier();
        if (marker_variance_.layout() == bayes::MarkerVarianceLayout::shared)
        {
            Eigen::Index n = 0;
            double sum_squares = 0.0;
            for (Eigen::Index i = 0; i < coeffs.size(); ++i)
            {
                const int component
                    = proportion == nullptr ? 1 : proportion->assignment(i);
                if (component == 0)
                {
                    continue;
                }
                ++n;
                if (multiplier == nullptr)
                {
                    sum_squares += coeffs(i) * coeffs(i);
                }
                else
                {
                    sum_squares
                        += (coeffs(i) * coeffs(i)) / (*multiplier)(component);
                }
            }
            if (n > 0)
            {
                variance_(0) = sampler({n, sum_squares}, rng_);
            }
            return;
        }

        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
        {
            const int component
                = proportion == nullptr ? 1 : proportion->assignment(i);
            if (component == 0)
            {
                continue;
            }
            double sum_squares = coeffs(i) * coeffs(i);
            if (multiplier != nullptr)
            {
                sum_squares /= (*multiplier)(component);
            }
            variance_(i) = sampler({1, sum_squares}, rng_);
        }
    }

   private:
    bayes::SingleGeneticBlockState& block_;
    const bayes::MarkerVariance& marker_variance_;
    Eigen::VectorXd& variance_;
    bayes::SingleProportionStateCap* proportion_cap_{};
    const bayes::SingleMultiplierCap* multiplier_cap_{};
    std::mt19937_64& rng_;
};

class JointGeneticVarianceStep final : public Step
{
   public:
    JointGeneticVarianceStep(
        const bayes::JointGeneticPrior& prior,
        bayes::JointGeneticBlockState& block,
        std::mt19937_64& rng)
        : prior_(prior), block_(block), rng_(rng)
    {
        if (prior_.query<bayes::JointMarkerVarianceCap>() == nullptr)
        {
            throw GelexException(
                "JointGeneticVarianceStep: prior lacks marker variance "
                "capability");
        }
        block_.prior_state().require<bayes::JointVarianceStateCap>();
    }

    auto step() -> void override
    {
        constexpr std::array modes{GeneticMode::A, GeneticMode::D};
        auto& variance_cap
            = block_.prior_state().require<bayes::JointVarianceStateCap>();
        auto* proportion_cap
            = block_.prior_state().query<bayes::JointProportionStateCap>();
        auto* proportion = proportion_cap == nullptr
                               ? nullptr
                               : &proportion_cap->proportion();
        for (auto mode : modes)
        {
            const auto& marker_variance
                = prior_.query<bayes::JointMarkerVarianceCap>()->variance(mode);
            auto& variance = variance_cap.variance(mode);
            const auto& coeffs = block_.state(mode).coeffs;
            auto sampler = make_variance_sampler(marker_variance.parameter());
            if (marker_variance.layout() == bayes::MarkerVarianceLayout::shared)
            {
                Eigen::Index n = 0;
                double sum_squares = 0.0;
                for (Eigen::Index i = 0; i < coeffs.size(); ++i)
                {
                    const int component
                        = proportion == nullptr ? 3 : proportion->assignment(i);
                    const bool active
                        = mode == GeneticMode::A
                              ? (component == 1 || component == 3)
                              : (component == 2 || component == 3);
                    if (!active)
                    {
                        continue;
                    }
                    ++n;
                    sum_squares += coeffs(i) * coeffs(i);
                }
                if (n > 0)
                {
                    variance(0) = sampler({n, sum_squares}, rng_);
                }
            }
            else
            {
                for (Eigen::Index i = 0; i < coeffs.size(); ++i)
                {
                    const int component
                        = proportion == nullptr ? 3 : proportion->assignment(i);
                    const bool active
                        = mode == GeneticMode::A
                              ? (component == 1 || component == 3)
                              : (component == 2 || component == 3);
                    if (active)
                    {
                        variance(i) = sampler({1, coeffs(i) * coeffs(i)}, rng_);
                    }
                }
            }
        }
    }

   private:
    const bayes::JointGeneticPrior& prior_;
    bayes::JointGeneticBlockState& block_;
    std::mt19937_64& rng_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEPS_GENETIC_VARIANCE_H_
