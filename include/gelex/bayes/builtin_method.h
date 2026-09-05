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

#ifndef GELEX_BAYES_BUILTIN_METHOD_H_
#define GELEX_BAYES_BUILTIN_METHOD_H_

#include <array>
#include <cstdint>
#include <string_view>
#include <utility>

#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/recipe.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

enum class BayesMethod : std::uint8_t
{
    RR,
    A,
    B,
    C,
    R,
    CD,
};

inline constexpr std::array bayes_method_names{
    std::pair{BayesMethod::RR, std::string_view{"RR"}},
    std::pair{BayesMethod::A, std::string_view{"A"}},
    std::pair{BayesMethod::B, std::string_view{"B"}},
    std::pair{BayesMethod::C, std::string_view{"C"}},
    std::pair{BayesMethod::R, std::string_view{"R"}},
    std::pair{BayesMethod::CD, std::string_view{"CD"}},
};

namespace detail
{

template <BayesMethod Method>
struct GeneticFamilyFor;

template <>
struct GeneticFamilyFor<BayesMethod::RR>
{
    using type = GaussianFamily<VarianceLayout::Pooled>;
};

template <>
struct GeneticFamilyFor<BayesMethod::A>
{
    using type = GaussianFamily<VarianceLayout::Unpooled>;
};

template <>
struct GeneticFamilyFor<BayesMethod::B>
{
    using type = SpikeSlabFamily<VarianceLayout::Unpooled>;
};

template <>
struct GeneticFamilyFor<BayesMethod::C>
{
    using type = SpikeSlabFamily<VarianceLayout::Pooled>;
};

template <>
struct GeneticFamilyFor<BayesMethod::R>
{
    using type = ScaledMixtureFamily<>;
};

template <>
struct GeneticFamilyFor<BayesMethod::CD>
{
    using type = JointSpikeSlabFamily<>;
};

}  // namespace detail

template <BayesMethod Method>
using genetic_family_t = typename detail::GeneticFamilyFor<Method>::type;

template <GeneticModeSet Modes, BayesMethod Method>
using BuiltinBayesRecipe = BayesRecipe<Modes, genetic_family_t<Method>>;

}  // namespace gelex

#endif  // GELEX_BAYES_BUILTIN_METHOD_H_
