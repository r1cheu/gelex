// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_BUILTIN_METHOD_H_
#define GELEX_BAYES_BUILTIN_METHOD_H_

#include <array>
#include <cstdint>
#include <string_view>
#include <utility>

#include "gelex/bayes/genetic/types.h"
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
    using type
        = HomogeneousModeValues<Modes, GaussianSpec<VarianceLayout::Pooled>>;
};

template <GeneticModeSet Modes>
struct BuiltinGeneticSpecFor<Modes, BayesMethod::A>
{
    using type
        = HomogeneousModeValues<Modes, GaussianSpec<VarianceLayout::Unpooled>>;
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
