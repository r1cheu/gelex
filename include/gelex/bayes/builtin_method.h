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

#include "gelex/bayes/genetic_policy.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
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

template <GeneticModeSet Modes, BayesMethod Method>
struct BuiltinGeneticSpecFor;

template <GeneticModeSet Modes>
struct BuiltinGeneticSpecFor<Modes, BayesMethod::RR>
{
    using type = GaussianSpec<VarianceLayout::Pooled>;
};

template <GeneticModeSet Modes>
struct BuiltinGeneticSpecFor<Modes, BayesMethod::A>
{
    using type = GaussianSpec<VarianceLayout::Unpooled>;
};

template <GeneticModeSet Modes>
struct BuiltinGeneticSpecFor<Modes, BayesMethod::B>
{
    using type
        = HomogeneousModeValues<Modes, SpikeSlabSpec<VarianceLayout::Unpooled>>;
};

template <GeneticModeSet Modes>
struct BuiltinGeneticSpecFor<Modes, BayesMethod::C>
{
    using type
        = HomogeneousModeValues<Modes, SpikeSlabSpec<VarianceLayout::Pooled>>;
};

template <GeneticModeSet Modes>
struct BuiltinGeneticSpecFor<Modes, BayesMethod::R>
{
    using type = HomogeneousModeValues<Modes, ScaledMixtureSpec<>>;
};

template <GeneticModeSet Modes>
    requires(Modes == (GeneticMode::A | GeneticMode::D))
struct BuiltinGeneticSpecFor<Modes, BayesMethod::CD>
{
    using type = JointModeValues<
        ModeValues<Modes, GaussianSpec<>, HalfNormalSpec>,
        JointSpikeSlabSpec<>>;
};

}  // namespace detail

template <GeneticModeSet Modes, BayesMethod Method>
using builtin_genetic_spec_t =
    typename detail::BuiltinGeneticSpecFor<Modes, Method>::type;

template <GeneticModeSet Modes, BayesMethod Method>
using BuiltinBayesRecipe
    = BayesRecipe<Modes, builtin_genetic_spec_t<Modes, Method>>;

}  // namespace gelex

#endif  // GELEX_BAYES_BUILTIN_METHOD_H_
