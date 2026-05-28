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

#ifndef GELEX_ALGO_INFER_MCMC_SWEEPS_GENETIC_COEFFICIENT_H_
#define GELEX_ALGO_INFER_MCMC_SWEEPS_GENETIC_COEFFICIENT_H_

#include <random>

#include <Eigen/Core>

#include "gelex/model/bayes/designs.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class SingleGeneticMarkerSweep
{
   public:
    explicit SingleGeneticMarkerSweep(const bayes::GeneticDesign& design)
        : design_(design)
    {
    }

    template <typename Transition>
    auto run(Transition& transition, std::mt19937_64& rng) -> void
    {
        const auto X = design_.X.matrix();
        const auto& XtX_diag = design_.XtX_diag;

        transition.prepare();
        for (Eigen::Index i = 0; i < X.cols(); ++i)
        {
            if (design_.is_monomorphic(i))
            {
                continue;
            }
            transition.update(i, X.col(i), XtX_diag(i), rng);
        }
        transition.finish();
    }

   private:
    const bayes::GeneticDesign& design_;
};

class JointGeneticMarkerSweep
{
   public:
    JointGeneticMarkerSweep(
        const bayes::GeneticDesign& additive,
        const bayes::GeneticDesign& dominance)
        : additive_(additive), dominance_(dominance)
    {
    }

    template <typename Transition>
    auto run(Transition& transition, std::mt19937_64& rng) -> void
    {
        const auto additive_x = additive_.X.matrix();
        const auto dominance_x = dominance_.X.matrix();

        transition.prepare();
        for (Eigen::Index i = 0; i < additive_x.cols(); ++i)
        {
            if (additive_.is_monomorphic(i) || dominance_.is_monomorphic(i))
            {
                continue;
            }
            transition.update(
                i,
                additive_x.col(i),
                additive_.XtX_diag(i),
                dominance_x.col(i),
                dominance_.XtX_diag(i),
                rng);
        }
        transition.finish();
    }

   private:
    const bayes::GeneticDesign& additive_;
    const bayes::GeneticDesign& dominance_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_SWEEPS_GENETIC_COEFFICIENT_H_
