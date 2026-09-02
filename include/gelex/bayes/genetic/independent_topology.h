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

#ifndef GELEX_BAYES_GENETIC_INDEPENDENT_TOPOLOGY_H_
#define GELEX_BAYES_GENETIC_INDEPENDENT_TOPOLOGY_H_

#include <array>
#include <concepts>
#include <cstddef>
#include <functional>
#include <ranges>
#include <type_traits>
#include <utility>

#include "gelex/types/genetic_mode.h"

namespace gelex
{

template <GeneticModeSet Modes, typename T>
class IndependentTopology
{
   public:
    using value_type = T;

    static constexpr GeneticModeSet modes = Modes;

    // Value-initializes every mode, so the defaults stay in T and are never
    // restated here. Implicitly deleted when T is not default constructible.
    constexpr IndependentTopology() = default;

    explicit constexpr IndependentTopology(T value)
        requires(Modes.size() == 1)
        : values_{std::move(value)}
    {
    }

    constexpr IndependentTopology(T additive, T dominance)
        requires(Modes == (GeneticMode::A | GeneticMode::D))
        : values_{std::move(additive), std::move(dominance)}
    {
    }

    // The one constructor that stays correct as modes are added; the named
    // overloads above exist only for readability at hand-written call sites.
    explicit constexpr IndependentTopology(std::array<T, Modes.size()> values)
        : values_(std::move(values))
    {
    }

    template <GeneticMode Mode>
    [[nodiscard]] constexpr auto get() noexcept -> T&
    {
        static constexpr std::size_t index = Modes.index_of(Mode);
        static_assert(index < Modes.size(), "mode not in topology");
        return values_[index];
    }

    template <GeneticMode Mode>
    [[nodiscard]] constexpr auto get() const noexcept -> const T&
    {
        static constexpr std::size_t index = Modes.index_of(Mode);
        static_assert(index < Modes.size(), "mode not in topology");
        return values_[index];
    }

    [[nodiscard]] constexpr auto each()
    {
        return std::views::zip(Modes.each(), values_);
    }

    [[nodiscard]] constexpr auto each() const
    {
        return std::views::zip(Modes.each(), values_);
    }

   private:
    std::array<T, Modes.size()> values_{};
};

template <GeneticModeSet Modes, typename Make>
    requires std::invocable<Make&, GeneticMode>
[[nodiscard]] constexpr auto generate_mode_values(Make make)
{
    using Value = std::remove_cvref_t<std::invoke_result_t<Make&, GeneticMode>>;
    return [&]<std::size_t... Index>(std::index_sequence<Index...>)
    {
        return IndependentTopology<Modes, Value>{
            std::array<Value, Modes.size()>{
                std::invoke(make, Modes.at(Index))...}};
    }(std::make_index_sequence<Modes.size()>{});
}

template <GeneticModeSet Modes, typename T, typename Convert>
    requires std::invocable<Convert&, GeneticMode, const T&>
[[nodiscard]] constexpr auto transform_mode_values(
    const IndependentTopology<Modes, T>& topology,
    Convert convert)
{
    using Result = std::remove_cvref_t<
        std::invoke_result_t<Convert&, GeneticMode, const T&>>;
    return [&]<std::size_t... Index>(std::index_sequence<Index...>)
    {
        return IndependentTopology<Modes, Result>{
            std::array<Result, Modes.size()>{std::invoke(
                convert,
                Modes.at(Index),
                topology.template get<Modes.at(Index)>())...}};
    }(std::make_index_sequence<Modes.size()>{});
}

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_INDEPENDENT_TOPOLOGY_H_
