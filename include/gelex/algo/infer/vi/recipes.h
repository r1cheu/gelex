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

#ifndef GELEX_ALGO_INFER_VI_RECIPES_H_
#define GELEX_ALGO_INFER_VI_RECIPES_H_

#include "gelex/algo/infer/chain.h"
#include "gelex/algo/infer/vi/context.h"
#include "gelex/algo/infer/vi/kernels/rr.h"
#include "gelex/algo/infer/vi/steps/fixed.h"
#include "gelex/algo/infer/vi/steps/genetic.h"
#include "gelex/algo/infer/vi/steps/random.h"
#include "gelex/algo/infer/vi/steps/residual.h"
#include "gelex/model/bayes/algorithm_shape.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::vi
{

template <bayes::AlgorithmShape Shape = bayes::AlgorithmShape::a_only>
inline auto make_bayes_rr_chain(const Context& ctx)
{
    static_assert(
        Shape != bayes::AlgorithmShape::ad_joint,
        "vi::make_bayes_rr_chain does not support ad_joint shape");

    if constexpr (Shape == bayes::AlgorithmShape::ad_independent)
    {
        return infer::Chain{
            FixedStep::make(ctx),
            RandomStep::make(ctx),
            GeneticStep<RRKernel>::make(ctx, GeneticMode::A),
            GeneticStep<RRKernel>::make(ctx, GeneticMode::D),
            ResidualStep::make(ctx),
        };
    }
    else
    {
        constexpr auto kMode = (Shape == bayes::AlgorithmShape::d_only)
                                   ? GeneticMode::D
                                   : GeneticMode::A;
        return infer::Chain{
            FixedStep::make(ctx),
            RandomStep::make(ctx),
            GeneticStep<RRKernel>::make(ctx, kMode),
            ResidualStep::make(ctx),
        };
    }
}

}  // namespace gelex::vi

#endif  // GELEX_ALGO_INFER_VI_RECIPES_H_
