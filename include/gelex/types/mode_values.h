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

#ifndef GELEX_TYPES_MODE_VALUES_H_
#define GELEX_TYPES_MODE_VALUES_H_

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

#include "gelex/types/genetic_mode.h"

namespace gelex
{

template <GeneticModeSet Modes, typename... Ts>
    requires(sizeof...(Ts) == Modes.size() && (std::is_object_v<Ts> && ...))
class ModeValues
{
    using values_type = std::tuple<Ts...>;

   public:
    static constexpr GeneticModeSet modes = Modes;

    template <GeneticMode Mode>
        requires(Modes.contains(Mode))
    using mode_value_type
        = std::tuple_element_t<Modes.index_of(Mode), values_type>;

    constexpr ModeValues() = default;

    explicit constexpr ModeValues(Ts... values) : values_{std::move(values)...}
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

template <GeneticMode Mode, typename Values, typename Convert>
constexpr auto transform_mode_value(const Values& values, Convert& convert)
    -> decltype(auto)
{
    return convert.template operator()<Mode>(values.template get<Mode>());
}

template <GeneticMode Mode, typename Values, typename Convert>
using transformed_mode_value_t
    = std::remove_cvref_t<decltype(transform_mode_value<Mode>(
        std::declval<const Values&>(),
        std::declval<Convert&>()))>;

}  // namespace detail

template <GeneticModeSet Modes, typename Make>
// NOLINTNEXTLINE(cppcoreguidelines-missing-std-forward)
[[nodiscard]] constexpr auto generate_mode_values(Make&& make)
{
    return [&]<std::size_t... Index>(std::index_sequence<Index...>)
    {
        return ModeValues<
            Modes,
            detail::generated_mode_value_t<Modes.at(Index), Make>...>{
            detail::generate_mode_value<Modes.at(Index)>(make)...};
    }(std::make_index_sequence<Modes.size()>{});
}

template <GeneticModeSet Modes, typename... Ts, typename Convert>
[[nodiscard]] constexpr auto transform_mode_values(
    const ModeValues<Modes, Ts...>& values,
    Convert&& convert)  // NOLINT(cppcoreguidelines-missing-std-forward)
{
    using values_type = ModeValues<Modes, Ts...>;
    return [&]<std::size_t... Index>(std::index_sequence<Index...>)
    {
        return ModeValues<
            Modes,
            detail::transformed_mode_value_t<
                Modes.at(Index),
                values_type,
                Convert>...>{
            detail::transform_mode_value<Modes.at(Index)>(values, convert)...};
    }(std::make_index_sequence<Modes.size()>{});
}

namespace detail
{

template <typename T, std::size_t>
using repeated_type_t = T;

template <GeneticModeSet Modes, typename T, std::size_t... Index>
auto homogeneous_mode_values_type(std::index_sequence<Index...>)
    -> ModeValues<Modes, repeated_type_t<T, Index>...>;

template <typename T>
inline constexpr bool is_mode_values = false;

template <GeneticModeSet Modes, typename... Ts>
inline constexpr bool is_mode_values<ModeValues<Modes, Ts...>> = true;

}  // namespace detail

template <GeneticModeSet Modes, typename T>
using HomogeneousModeValues
    = decltype(detail::homogeneous_mode_values_type<Modes, T>(
        std::make_index_sequence<Modes.size()>{}));

template <typename ModeValuesType, typename JointT>
    requires(
        detail::is_mode_values<ModeValuesType>
        && ModeValuesType::modes == (GeneticMode::A | GeneticMode::D)
        && std::is_object_v<JointT>)
class JointModeValues
{
   public:
    using mode_values_type = ModeValuesType;
    using joint_value_type = JointT;

    template <GeneticMode Mode>
    using mode_value_type =
        typename mode_values_type::template mode_value_type<Mode>;

    static constexpr GeneticModeSet modes = mode_values_type::modes;

    constexpr JointModeValues() = default;

    constexpr JointModeValues(ModeValuesType mode_values, JointT joint)
        : mode_values_{std::move(mode_values)}, joint_{std::move(joint)}
    {
    }

    [[nodiscard]] constexpr auto mode_values() noexcept -> mode_values_type&
    {
        return mode_values_;
    }

    [[nodiscard]] constexpr auto mode_values() const noexcept
        -> const mode_values_type&
    {
        return mode_values_;
    }

    [[nodiscard]] constexpr auto joint() noexcept -> JointT& { return joint_; }

    [[nodiscard]] constexpr auto joint() const noexcept -> const JointT&
    {
        return joint_;
    }

   private:
    mode_values_type mode_values_;
    JointT joint_;
};

template <typename ModeValuesType, typename JointT>
JointModeValues(ModeValuesType, JointT)
    -> JointModeValues<ModeValuesType, JointT>;

}  // namespace gelex

#endif  // GELEX_TYPES_MODE_VALUES_H_
