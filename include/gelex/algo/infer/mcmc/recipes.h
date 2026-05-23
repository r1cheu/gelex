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

#ifndef GELEX_ALGO_INFER_MCMC_RECIPES_H_
#define GELEX_ALGO_INFER_MCMC_RECIPES_H_

#include "gelex/algo/infer/chain.h"
#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_a.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_b.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_c.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_r.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_rr.h"
#include "gelex/algo/infer/mcmc/steps/fixed.h"
#include "gelex/algo/infer/mcmc/steps/genetic.h"
#include "gelex/algo/infer/mcmc/steps/genetic_joint.h"
#include "gelex/algo/infer/mcmc/steps/pi.h"
#include "gelex/algo/infer/mcmc/steps/random.h"
#include "gelex/algo/infer/mcmc/steps/residual.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

enum class GeneticShape
{
    a_only,
    d_only,
    ad_independent,
};

template <typename Kernel, GeneticShape Shape = GeneticShape::a_only>
inline auto make_simple_chain(const Context& ctx)
{
    if constexpr (Shape == GeneticShape::ad_independent)
    {
        return infer::Chain{
            FixedStep::make(ctx),
            RandomStep::make(ctx),
            GeneticStep<Kernel>::make(ctx, GeneticMode::A),
            GeneticStep<Kernel>::make(ctx, GeneticMode::D),
            ResidualStep::make(ctx),
        };
    }
    else
    {
        constexpr auto kMode
            = (Shape == GeneticShape::d_only) ? GeneticMode::D : GeneticMode::A;
        return infer::Chain{
            FixedStep::make(ctx),
            RandomStep::make(ctx),
            GeneticStep<Kernel>::make(ctx, kMode),
            ResidualStep::make(ctx),
        };
    }
}

template <typename Kernel, GeneticShape Shape = GeneticShape::a_only>
inline auto make_pi_chain(const Context& ctx)
{
    if constexpr (Shape == GeneticShape::ad_independent)
    {
        return infer::Chain{
            FixedStep::make(ctx),
            RandomStep::make(ctx),
            GeneticStep<Kernel>::make(ctx, GeneticMode::A),
            PiStep::make(ctx, GeneticMode::A),
            GeneticStep<Kernel>::make(ctx, GeneticMode::D),
            PiStep::make(ctx, GeneticMode::D),
            ResidualStep::make(ctx),
        };
    }
    else
    {
        constexpr auto kMode
            = (Shape == GeneticShape::d_only) ? GeneticMode::D : GeneticMode::A;
        return infer::Chain{
            FixedStep::make(ctx),
            RandomStep::make(ctx),
            GeneticStep<Kernel>::make(ctx, kMode),
            PiStep::make(ctx, kMode),
            ResidualStep::make(ctx),
        };
    }
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_a_chain(const Context& ctx)
{
    return make_simple_chain<BayesAKernel, Shape>(ctx);
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_b_chain(const Context& ctx)
{
    return make_simple_chain<BayesBKernel, Shape>(ctx);
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_bpi_chain(const Context& ctx)
{
    return make_pi_chain<BayesBKernel, Shape>(ctx);
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_c_chain(const Context& ctx)
{
    return make_simple_chain<BayesCKernel, Shape>(ctx);
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_cpi_chain(const Context& ctx)
{
    return make_pi_chain<BayesCKernel, Shape>(ctx);
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_r_chain(const Context& ctx)
{
    return make_pi_chain<BayesRKernel, Shape>(ctx);
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_rr_chain(const Context& ctx)
{
    return make_simple_chain<BayesRRKernel, Shape>(ctx);
}

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_RECIPES_H_
