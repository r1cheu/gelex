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

#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/exception.h"

namespace gelex::mcmc
{

enum class GeneticShape
{
    a_only,
    d_only,
    ad_independent,
};

class UnsupportedChain
{
   public:
    auto step() -> void
    {
        throw GelexException(
            "MCMC recipe chain is not implemented after Bayes prior/state "
            "cleanup");
    }
};

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_a_chain(const Context& ctx)
{
    static_cast<void>(ctx);
    static_cast<void>(Shape);
    return UnsupportedChain{};
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_b_chain(const Context& ctx)
{
    static_cast<void>(ctx);
    static_cast<void>(Shape);
    return UnsupportedChain{};
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_bpi_chain(const Context& ctx)
{
    static_cast<void>(ctx);
    static_cast<void>(Shape);
    return UnsupportedChain{};
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_c_chain(const Context& ctx)
{
    static_cast<void>(ctx);
    static_cast<void>(Shape);
    return UnsupportedChain{};
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_cpi_chain(const Context& ctx)
{
    static_cast<void>(ctx);
    static_cast<void>(Shape);
    return UnsupportedChain{};
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_r_chain(const Context& ctx)
{
    static_cast<void>(ctx);
    static_cast<void>(Shape);
    return UnsupportedChain{};
}

template <GeneticShape Shape = GeneticShape::a_only>
inline auto make_bayes_rr_chain(const Context& ctx)
{
    static_cast<void>(ctx);
    static_cast<void>(Shape);
    return UnsupportedChain{};
}

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_RECIPES_H_
