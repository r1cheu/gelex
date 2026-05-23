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

#ifndef GELEX_MODEL_BAYES_GENETIC_PRIOR_H_
#define GELEX_MODEL_BAYES_GENETIC_PRIOR_H_

#include <algorithm>
#include <memory>
#include <span>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/prior_specs.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

class GeneticPrior
{
   public:
    auto operator=(const GeneticPrior&) -> GeneticPrior& = delete;
    auto operator=(GeneticPrior&&) noexcept -> GeneticPrior& = delete;

    virtual ~GeneticPrior() = default;

    virtual auto modes() const -> std::span<const GeneticMode> = 0;

    virtual auto make_state(
        Eigen::Index num_markers,
        Eigen::Index num_individuals) const
        -> std::unique_ptr<GeneticPriorState>
        = 0;

    auto contains(GeneticMode mode) const -> bool
    {
        return std::ranges::contains(modes(), mode);
    }

    template <typename Capability>
    auto query() const -> const Capability*
    {
        return dynamic_cast<const Capability*>(this);
    }

   protected:
    GeneticPrior() = default;
    GeneticPrior(const GeneticPrior&) = default;
    GeneticPrior(GeneticPrior&&) noexcept = default;

    static auto make_variance_values(
        std::span<const MarkerVarianceSpec> variance_specs,
        Eigen::Index num_markers) -> std::vector<Eigen::VectorXd>;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GENETIC_PRIOR_H_
