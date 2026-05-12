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

#ifndef GELEX_MODEL_BAYES_GENETIC_PRIOR_RUNTIME_STATE_H_
#define GELEX_MODEL_BAYES_GENETIC_PRIOR_RUNTIME_STATE_H_

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

struct GeneticEffectRuntimeInit
{
    GeneticMode mode{};
    Eigen::Index num_markers{};
};

// effects[i].mode == prior.modes()[i]; built by BayesState::init() after
// validate_prior_for_model() — concrete priors index, not search.
struct GeneticPriorRuntimeInit
{
    std::span<const GeneticEffectRuntimeInit> effects;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GENETIC_PRIOR_RUNTIME_STATE_H_
