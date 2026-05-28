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
#include <type_traits>
#include <variant>

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
class SingleSharedGeneticVarianceStep final : public Step
{
   public:
    SingleSharedGeneticVarianceStep(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        std::mt19937_64& rng)
        : block_(block),
          variance_parameter_(
              prior.query<bayes::SingleSharedMarkerVarianceCap>() != nullptr
                  ? prior.query<bayes::SingleSharedMarkerVarianceCap>()
                        ->variance()
                        .parameter()
                  : throw GelexException(
                        "SingleSharedGeneticVarianceStep: prior lacks shared "
                        "variance capability")),
          variance_(block.prior_state()
                        .require<bayes::SingleSharedVarianceStateCap>()
                        .variance()),
          proportion_cap_(
              block.prior_state().query<bayes::SingleProportionStateCap>()),
          multiplier_cap_(prior.query<bayes::SingleMultiplierCap>()),
          rng_(rng)
    {
    }

    auto step() -> void override
    {
        auto sampler = make_variance_sampler(variance_parameter_);
        const auto& coeffs = block_.state().coeffs;
        auto* proportion = proportion_cap_ == nullptr
                               ? nullptr
                               : &proportion_cap_->proportion();
        const auto* multiplier = multiplier_cap_ == nullptr
                                     ? nullptr
                                     : &multiplier_cap_->multiplier();
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
            variance_ = sampler({n, sum_squares}, rng_);
        }
    }

   private:
    bayes::SingleGeneticBlockState& block_;
    const bayes::VarianceParameter& variance_parameter_;
    double& variance_;
    bayes::SingleProportionStateCap* proportion_cap_{};
    const bayes::SingleMultiplierCap* multiplier_cap_{};
    std::mt19937_64& rng_;
};

class SinglePerMarkerGeneticVarianceStep final : public Step
{
   public:
    SinglePerMarkerGeneticVarianceStep(
        const bayes::SingleGeneticPrior& prior,
        bayes::SingleGeneticBlockState& block,
        std::mt19937_64& rng)
        : block_(block),
          variance_parameter_(
              prior.query<bayes::SinglePerMarkerVarianceCap>() != nullptr
                  ? prior.query<bayes::SinglePerMarkerVarianceCap>()
                        ->variance()
                        .parameter()
                  : throw GelexException(
                        "SinglePerMarkerGeneticVarianceStep: prior lacks "
                        "per-marker variance capability")),
          variance_(block.prior_state()
                        .require<bayes::SinglePerMarkerVarianceStateCap>()
                        .variance()),
          proportion_cap_(
              block.prior_state().query<bayes::SingleProportionStateCap>()),
          rng_(rng)
    {
    }

    auto step() -> void override
    {
        auto sampler = make_variance_sampler(variance_parameter_);
        const auto& coeffs = block_.state().coeffs;
        auto* proportion = proportion_cap_ == nullptr
                               ? nullptr
                               : &proportion_cap_->proportion();
        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
        {
            const int component
                = proportion == nullptr ? 1 : proportion->assignment(i);
            if (component == 0)
            {
                continue;
            }
            double sum_squares = coeffs(i) * coeffs(i);
            variance_(i) = sampler({1, sum_squares}, rng_);
        }
    }

   private:
    bayes::SingleGeneticBlockState& block_;
    const bayes::VarianceParameter& variance_parameter_;
    Eigen::VectorXd& variance_;
    bayes::SingleProportionStateCap* proportion_cap_{};
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
            std::visit(
                [&](const auto& marker_variance_slot)
                {
                    auto sampler = make_variance_sampler(
                        marker_variance_slot.parameter());
                    using MarkerVarianceType
                        = std::decay_t<decltype(marker_variance_slot)>;
                    if constexpr (
                        std::is_same_v<
                            MarkerVarianceType,
                            bayes::SharedMarkerVariance>)
                    {
                        auto& shared_variance = std::get<double>(variance);
                        Eigen::Index n = 0;
                        double sum_squares = 0.0;
                        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
                        {
                            const int component
                                = proportion == nullptr
                                      ? 3
                                      : proportion->assignment(i);
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
                            shared_variance = sampler({n, sum_squares}, rng_);
                        }
                    }
                    else
                    {
                        auto& per_marker_variance
                            = std::get<Eigen::VectorXd>(variance);
                        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
                        {
                            const int component
                                = proportion == nullptr
                                      ? 3
                                      : proportion->assignment(i);
                            const bool active
                                = mode == GeneticMode::A
                                      ? (component == 1 || component == 3)
                                      : (component == 2 || component == 3);
                            if (active)
                            {
                                per_marker_variance(i)
                                    = sampler({1, coeffs(i) * coeffs(i)}, rng_);
                            }
                        }
                    }
                },
                marker_variance);
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
