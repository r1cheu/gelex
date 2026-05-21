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

#ifndef GELEX_MODEL_BAYES_RUNTIME_STATE_H_
#define GELEX_MODEL_BAYES_RUNTIME_STATE_H_

#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

class GeneticPriorRuntimeState
{
   public:
    auto operator=(const GeneticPriorRuntimeState&)
        -> GeneticPriorRuntimeState& = delete;
    auto operator=(GeneticPriorRuntimeState&&) noexcept
        -> GeneticPriorRuntimeState& = delete;

    virtual ~GeneticPriorRuntimeState() = default;

   protected:
    GeneticPriorRuntimeState() = default;
    GeneticPriorRuntimeState(const GeneticPriorRuntimeState&) = default;
    GeneticPriorRuntimeState(GeneticPriorRuntimeState&&) noexcept = default;
};

class MarkerVarianceRuntimeState
{
   public:
    explicit MarkerVarianceRuntimeState(std::vector<Eigen::VectorXd> variances)
        : variances_(std::move(variances))
    {
    }

    auto variances() -> std::span<Eigen::VectorXd> { return variances_; }
    auto variances() const -> std::span<const Eigen::VectorXd>
    {
        return variances_;
    }

    auto variance_at(std::size_t i) -> Eigen::VectorXd&
    {
        return variances_[i];
    }
    auto variance_at(std::size_t i) const -> const Eigen::VectorXd&
    {
        return variances_[i];
    }

   private:
    std::vector<Eigen::VectorXd> variances_;
};

class MixtureRuntimeState
{
   public:
    MixtureRuntimeState(Eigen::VectorXi assignment, Eigen::VectorXd proportion)
        : assignment_(std::move(assignment)), proportion_(std::move(proportion))
    {
    }

    auto assignment() -> Eigen::VectorXi& { return assignment_; }
    auto assignment() const -> const Eigen::VectorXi& { return assignment_; }

    auto proportion() -> Eigen::VectorXd& { return proportion_; }
    auto proportion() const -> const Eigen::VectorXd& { return proportion_; }

   private:
    Eigen::VectorXi assignment_;
    Eigen::VectorXd proportion_;
};

class GaussianRuntimeState final : public GeneticPriorRuntimeState
{
   public:
    explicit GaussianRuntimeState(MarkerVarianceRuntimeState marker_variance)
        : marker_variance_(std::move(marker_variance))
    {
    }

    auto marker_variance() -> MarkerVarianceRuntimeState&
    {
        return marker_variance_;
    }
    auto marker_variance() const -> const MarkerVarianceRuntimeState&
    {
        return marker_variance_;
    }

   private:
    MarkerVarianceRuntimeState marker_variance_;
};

class MixtureGaussianRuntimeState final : public GeneticPriorRuntimeState
{
   public:
    MixtureGaussianRuntimeState(
        MarkerVarianceRuntimeState marker_variance,
        MixtureRuntimeState mixture)
        : marker_variance_(std::move(marker_variance)),
          mixture_(std::move(mixture))
    {
    }

    auto marker_variance() -> MarkerVarianceRuntimeState&
    {
        return marker_variance_;
    }
    auto marker_variance() const -> const MarkerVarianceRuntimeState&
    {
        return marker_variance_;
    }

    auto mixture() -> MixtureRuntimeState& { return mixture_; }
    auto mixture() const -> const MixtureRuntimeState& { return mixture_; }

   private:
    MarkerVarianceRuntimeState marker_variance_;
    MixtureRuntimeState mixture_;
};

struct GeneticEffectRuntimeInit
{
    GeneticMode mode{};
    Eigen::Index num_markers{};
};

struct GeneticPriorRuntimeInit
{
    std::span<const GeneticEffectRuntimeInit> effects;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_RUNTIME_STATE_H_
