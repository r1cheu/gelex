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

#include "gelex/bayes/semantic_method.h"

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

inline constexpr std::array BAYES_METHOD_NAMES{
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
struct SemanticMethodFor;

template <>
struct SemanticMethodFor<BayesMethod::RR>
{
    using type = GaussianMethod<Variance::Pooled>;
};

template <>
struct SemanticMethodFor<BayesMethod::A>
{
    using type = GaussianMethod<Variance::Unpooled>;
};

template <>
struct SemanticMethodFor<BayesMethod::B>
{
    using type = SpikeSlabMethod<Variance::Unpooled>;
};

template <>
struct SemanticMethodFor<BayesMethod::C>
{
    using type = SpikeSlabMethod<Variance::Pooled>;
};

template <>
struct SemanticMethodFor<BayesMethod::R>
{
    using type = ScaledMixtureMethod;
};

template <>
struct SemanticMethodFor<BayesMethod::CD>
{
    using type = JointSpikeSlabMethod;
};

}  // namespace detail

template <BayesMethod Method>
using semantic_method_t = typename detail::SemanticMethodFor<Method>::type;

}  // namespace gelex

#endif  // GELEX_BAYES_BUILTIN_METHOD_H_
