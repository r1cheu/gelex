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

#include "gelex/algo/infer/mcmc/chain.h"
#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_a.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_b.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_c.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_r.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_rr.h"
#include "gelex/algo/infer/mcmc/samplers/fixed.h"
#include "gelex/algo/infer/mcmc/samplers/genetic.h"
#include "gelex/algo/infer/mcmc/samplers/pi.h"
#include "gelex/algo/infer/mcmc/samplers/random.h"
#include "gelex/algo/infer/mcmc/samplers/residual.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

inline auto make_bayes_a_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesAKernel>::make(ctx, GeneticMode::A),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_ad_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesAKernel>::make(ctx, GeneticMode::A),
        GeneticSampler<BayesAKernel>::make(ctx, GeneticMode::D),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_b_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesBKernel>::make(ctx, GeneticMode::A),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_bpi_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesBKernel>::make(ctx, GeneticMode::A),
        PiSampler::make(ctx, GeneticMode::A),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_bd_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesBKernel>::make(ctx, GeneticMode::A),
        GeneticSampler<BayesBKernel>::make(ctx, GeneticMode::D),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_bdpi_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesBKernel>::make(ctx, GeneticMode::A),
        PiSampler::make(ctx, GeneticMode::A),
        GeneticSampler<BayesBKernel>::make(ctx, GeneticMode::D),
        PiSampler::make(ctx, GeneticMode::D),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_c_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesCKernel>::make(ctx, GeneticMode::A),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_cpi_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesCKernel>::make(ctx, GeneticMode::A),
        PiSampler::make(ctx, GeneticMode::A),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_cd_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesCKernel>::make(ctx, GeneticMode::A),
        GeneticSampler<BayesCKernel>::make(ctx, GeneticMode::D),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_cdpi_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesCKernel>::make(ctx, GeneticMode::A),
        PiSampler::make(ctx, GeneticMode::A),
        GeneticSampler<BayesCKernel>::make(ctx, GeneticMode::D),
        PiSampler::make(ctx, GeneticMode::D),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_rr_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesRRKernel>::make(ctx, GeneticMode::A),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_rrd_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesRRKernel>::make(ctx, GeneticMode::A),
        GeneticSampler<BayesRRKernel>::make(ctx, GeneticMode::D),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_r_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesRKernel>::make(ctx, GeneticMode::A),
        PiSampler::make(ctx, GeneticMode::A),
        ResidualSampler::make(ctx),
    };
}

inline auto make_bayes_rd_chain(const Context& ctx)
{
    return McmcChain{
        FixedSampler::make(ctx),
        RandomSampler::make(ctx),
        GeneticSampler<BayesRKernel>::make(ctx, GeneticMode::A),
        PiSampler::make(ctx, GeneticMode::A),
        GeneticSampler<BayesRKernel>::make(ctx, GeneticMode::D),
        PiSampler::make(ctx, GeneticMode::D),
        ResidualSampler::make(ctx),
    };
}

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_RECIPES_H_
