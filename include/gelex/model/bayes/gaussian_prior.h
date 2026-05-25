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

#ifndef GELEX_MODEL_BAYES_GAUSSIAN_PRIOR_H_
#define GELEX_MODEL_BAYES_GAUSSIAN_PRIOR_H_

#include <array>
#include <memory>
#include <span>
#include <string_view>

#include <Eigen/Core>

#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

class GaussianPrior final
    : public GeneticPrior
    , public MarkerVarianceCap
{
   public:
    static constexpr std::string_view name = "gaussian";

    GaussianPrior(GeneticMode mode, MarkerVariance variance);

    auto modes() const -> std::span<const GeneticMode> override
    {
        return modes_;
    }
    auto variance() -> std::span<MarkerVariance> override
    {
        return marker_variances_;
    }
    auto variance() const -> std::span<const MarkerVariance> override
    {
        return marker_variances_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<GeneticPriorState> override;

   private:
    std::array<GeneticMode, 1> modes_;
    std::array<MarkerVariance, 1> marker_variances_;
};

class SpikeSlabGaussianPrior final
    : public GeneticPrior
    , public MarkerVarianceCap
    , public MixtureProportionCap
{
   public:
    static constexpr std::string_view name = "spike_slab_gaussian";

    SpikeSlabGaussianPrior(
        GeneticMode mode,
        MarkerVariance variance,
        MixtureProportion proportion);

    auto modes() const -> std::span<const GeneticMode> override
    {
        return modes_;
    }
    auto variance() -> std::span<MarkerVariance> override
    {
        return marker_variances_;
    }
    auto variance() const -> std::span<const MarkerVariance> override
    {
        return marker_variances_;
    }
    auto proportion() -> std::span<MixtureProportion> override
    {
        return mixture_proportions_;
    }
    auto proportion() const -> std::span<const MixtureProportion> override
    {
        return mixture_proportions_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<GeneticPriorState> override;

   private:
    std::array<GeneticMode, 1> modes_;
    std::array<MarkerVariance, 1> marker_variances_;
    std::array<MixtureProportion, 1> mixture_proportions_;
};

class ScaledMixtureGaussianPrior final
    : public GeneticPrior
    , public MarkerVarianceCap
    , public MixtureProportionCap
    , public MultiplierCap
{
   public:
    static constexpr std::string_view name = "scaled_mixture_gaussian";

    ScaledMixtureGaussianPrior(
        GeneticMode mode,
        MarkerVariance variance,
        Eigen::VectorXd multiplier,
        MixtureProportion proportion);

    auto modes() const -> std::span<const GeneticMode> override
    {
        return modes_;
    }
    auto variance() -> std::span<MarkerVariance> override
    {
        return marker_variances_;
    }
    auto variance() const -> std::span<const MarkerVariance> override
    {
        return marker_variances_;
    }
    auto proportion() -> std::span<MixtureProportion> override
    {
        return mixture_proportions_;
    }
    auto proportion() const -> std::span<const MixtureProportion> override
    {
        return mixture_proportions_;
    }
    auto multiplier() -> std::span<Eigen::VectorXd> override
    {
        return multipliers_;
    }
    auto multiplier() const -> std::span<const Eigen::VectorXd> override
    {
        return multipliers_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<GeneticPriorState> override;

   private:
    std::array<GeneticMode, 1> modes_;
    std::array<MarkerVariance, 1> marker_variances_;
    std::array<Eigen::VectorXd, 1> multipliers_;
    std::array<MixtureProportion, 1> mixture_proportions_;
};

class JointMixtureGaussianPrior final
    : public GeneticPrior
    , public MarkerVarianceCap
    , public MixtureProportionCap
{
   public:
    static constexpr std::string_view name = "joint_mixture_gaussian";

    JointMixtureGaussianPrior(
        std::array<GeneticMode, 2> modes,
        std::array<MarkerVariance, 2> variances,
        MixtureProportion proportion);

    auto modes() const -> std::span<const GeneticMode> override
    {
        return modes_;
    }
    auto variance() -> std::span<MarkerVariance> override
    {
        return marker_variances_;
    }
    auto variance() const -> std::span<const MarkerVariance> override
    {
        return marker_variances_;
    }
    auto proportion() -> std::span<MixtureProportion> override
    {
        return mixture_proportions_;
    }
    auto proportion() const -> std::span<const MixtureProportion> override
    {
        return mixture_proportions_;
    }

    auto visit(infra::FieldVisitor& visitor) -> void override;
    auto make_state(Eigen::Index num_markers, Eigen::Index num_individuals)
        const -> std::unique_ptr<GeneticPriorState> override;

   private:
    std::array<GeneticMode, 2> modes_;
    std::array<MarkerVariance, 2> marker_variances_;
    std::array<MixtureProportion, 1> mixture_proportions_;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GAUSSIAN_PRIOR_H_
