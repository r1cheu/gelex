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

#include <Eigen/Core>

#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior_capabilities.h"
#include "gelex/model/bayes/prior_specs.h"
#include "gelex/model/bayes/runtime_state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

class GaussianPrior final
    : public GeneticPrior
    , public VarianceCapability
{
   public:
    GaussianPrior(GeneticMode mode, MarkerVarianceSpec variance);

    auto modes() const -> std::span<const GeneticMode> override;
    auto variance_specs() const -> std::span<const MarkerVarianceSpec> override;
    auto make_state(const GeneticPriorRuntimeInit& init) const
        -> std::unique_ptr<GeneticPriorRuntimeState> override;

   private:
    std::array<GeneticMode, 1> modes_;
    std::array<MarkerVarianceSpec, 1> variance_specs_;
};

class SpikeSlabGaussianPrior final
    : public GeneticPrior
    , public VarianceCapability
    , public ProportionCapability
{
   public:
    SpikeSlabGaussianPrior(
        GeneticMode mode,
        MarkerVarianceSpec variance,
        ProportionSpec proportion);

    auto modes() const -> std::span<const GeneticMode> override;
    auto variance_specs() const -> std::span<const MarkerVarianceSpec> override;
    auto proportion_specs() const -> std::span<const ProportionSpec> override;
    auto make_state(const GeneticPriorRuntimeInit& init) const
        -> std::unique_ptr<GeneticPriorRuntimeState> override;

   private:
    std::array<GeneticMode, 1> modes_;
    std::array<MarkerVarianceSpec, 1> variance_specs_;
    std::array<ProportionSpec, 1> proportion_specs_;
};

class ScaledMixtureGaussianPrior final
    : public GeneticPrior
    , public VarianceCapability
    , public ProportionCapability
    , public MultiplierCapability
{
   public:
    ScaledMixtureGaussianPrior(
        GeneticMode mode,
        MarkerVarianceSpec variance,
        Eigen::VectorXd multiplier,
        ProportionSpec proportion);

    auto modes() const -> std::span<const GeneticMode> override;
    auto variance_specs() const -> std::span<const MarkerVarianceSpec> override;
    auto proportion_specs() const -> std::span<const ProportionSpec> override;
    auto multipliers() const -> std::span<const Eigen::VectorXd> override;
    auto make_state(const GeneticPriorRuntimeInit& init) const
        -> std::unique_ptr<GeneticPriorRuntimeState> override;

   private:
    std::array<GeneticMode, 1> modes_;
    std::array<MarkerVarianceSpec, 1> variance_specs_;
    std::array<Eigen::VectorXd, 1> multipliers_;
    std::array<ProportionSpec, 1> proportion_specs_;
};

class JointMixtureGaussianPrior final
    : public GeneticPrior
    , public VarianceCapability
    , public ProportionCapability
{
   public:
    JointMixtureGaussianPrior(
        std::array<GeneticMode, 2> modes,
        std::array<MarkerVarianceSpec, 2> variances,
        ProportionSpec proportion);

    auto modes() const -> std::span<const GeneticMode> override;
    auto variance_specs() const -> std::span<const MarkerVarianceSpec> override;
    auto proportion_specs() const -> std::span<const ProportionSpec> override;
    auto make_state(const GeneticPriorRuntimeInit& init) const
        -> std::unique_ptr<GeneticPriorRuntimeState> override;

   private:
    std::array<GeneticMode, 2> modes_;
    std::array<MarkerVarianceSpec, 2> variance_specs_;
    std::array<ProportionSpec, 1> proportion_specs_;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GAUSSIAN_PRIOR_H_
