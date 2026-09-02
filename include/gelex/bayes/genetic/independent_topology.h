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

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

#include "gelex/types/genetic_mode.h"

namespace gelex
{

template <GeneticModeSet Modes, typename... Ts>
    requires(sizeof...(Ts) == Modes.size() && (std::is_object_v<Ts> && ...))
class IndependentTopology
{
    using values_type = std::tuple<Ts...>;

   public:
    static constexpr GeneticModeSet modes = Modes;

    template <GeneticMode Mode>
        requires(Modes.contains(Mode))
    using mode_value_type
        = std::tuple_element_t<Modes.index_of(Mode), values_type>;

    constexpr IndependentTopology() = default;

    explicit constexpr IndependentTopology(Ts... values)
        : values_{std::move(values)...}
    {
    }

    template <GeneticMode Mode>
        requires(Modes.contains(Mode))
    [[nodiscard]] constexpr auto get() noexcept -> mode_value_type<Mode>&
    {
        return std::get<Modes.index_of(Mode)>(values_);
    }

    template <GeneticMode Mode>
        requires(Modes.contains(Mode))
    [[nodiscard]] constexpr auto get() const noexcept
        -> const mode_value_type<Mode>&
    {
        return std::get<Modes.index_of(Mode)>(values_);
    }

    template <typename Function>
    // NOLINTNEXTLINE(cppcoreguidelines-missing-std-forward)
    constexpr auto for_each(Function&& function) -> void
    {
        [&]<std::size_t... Index>(std::index_sequence<Index...>)
        {
            (function.template operator()<Modes.at(Index)>(
                 this->template get<Modes.at(Index)>()),
             ...);
        }(std::make_index_sequence<Modes.size()>{});
    }

    template <typename Function>
    // NOLINTNEXTLINE(cppcoreguidelines-missing-std-forward)
    constexpr auto for_each(Function&& function) const -> void
    {
        [&]<std::size_t... Index>(std::index_sequence<Index...>)
        {
            (function.template operator()<Modes.at(Index)>(
                 this->template get<Modes.at(Index)>()),
             ...);
        }(std::make_index_sequence<Modes.size()>{});
    }

   private:
    std::tuple<Ts...> values_;
};

namespace detail
{

template <GeneticMode Mode, typename Make>
constexpr auto generate_mode_value(Make& make) -> decltype(auto)
{
    return make.template operator()<Mode>();
}

template <GeneticMode Mode, typename Make>
using generated_mode_value_t
    = std::remove_cvref_t<decltype(generate_mode_value<Mode>(
        std::declval<Make&>()))>;

template <GeneticMode Mode, typename Topology, typename Convert>
constexpr auto transform_mode_value(const Topology& topology, Convert& convert)
    -> decltype(auto)
{
    return convert.template operator()<Mode>(topology.template get<Mode>());
}

template <GeneticMode Mode, typename Topology, typename Convert>
using transformed_mode_value_t
    = std::remove_cvref_t<decltype(transform_mode_value<Mode>(
        std::declval<const Topology&>(),
        std::declval<Convert&>()))>;

}  // namespace detail

template <GeneticModeSet Modes, typename Make>
// NOLINTNEXTLINE(cppcoreguidelines-missing-std-forward)
[[nodiscard]] constexpr auto generate_mode_values(Make&& make)
{
    return [&]<std::size_t... Index>(std::index_sequence<Index...>)
    {
        return IndependentTopology<
            Modes,
            detail::generated_mode_value_t<Modes.at(Index), Make>...>{
            detail::generate_mode_value<Modes.at(Index)>(make)...};
    }(std::make_index_sequence<Modes.size()>{});
}

template <GeneticModeSet Modes, typename... Ts, typename Convert>
[[nodiscard]] constexpr auto transform_mode_values(
    const IndependentTopology<Modes, Ts...>& topology,
    Convert&& convert)  // NOLINT(cppcoreguidelines-missing-std-forward)
{
    using topology_type = IndependentTopology<Modes, Ts...>;
    return [&]<std::size_t... Index>(std::index_sequence<Index...>)
    {
        return IndependentTopology<
            Modes,
            detail::transformed_mode_value_t<
                Modes.at(Index),
                topology_type,
                Convert>...>{
            detail::transform_mode_value<Modes.at(Index)>(
                topology, convert)...};
    }(std::make_index_sequence<Modes.size()>{});
}

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_INDEPENDENT_TOPOLOGY_H_
