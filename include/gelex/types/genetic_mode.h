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

#ifndef GELEX_TYPES_GENETIC_MODE_H_
#define GELEX_TYPES_GENETIC_MODE_H_
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <iterator>
#include <limits>
#include <map>
#include <ranges>
#include <string_view>
#include <utility>

namespace gelex
{

enum class GeneticMode : std::uint8_t
{
    A,
    D,
};

inline constexpr std::array genetic_mode_names{
    std::pair{GeneticMode::A, std::string_view{"A"}},
    std::pair{GeneticMode::D, std::string_view{"D"}},
};

inline constexpr std::array all_genetic_modes{GeneticMode::A, GeneticMode::D};

template <typename T>
using ModeMap = std::map<GeneticMode, T>;

// A non-empty subset of GeneticMode. each() yields its members in enum order.
class GeneticModeSet
{
   public:
    // Implementation detail exposed only because a non-type template parameter
    // requires a structural type; use the member functions instead. Undefined
    // bits are ignored by every observer, so a stray write cannot desync them.
    std::uint8_t bits{};

    explicit constexpr GeneticModeSet(GeneticMode mode) noexcept
        : bits(bit_for(mode))
    {
    }

    [[nodiscard]] constexpr auto contains(GeneticMode mode) const noexcept
        -> bool
    {
        return (bits & bit_for(mode)) != 0;
    }

    [[nodiscard]] constexpr auto each() const
    {
        return all_genetic_modes
               | std::views::filter([bits = this->bits](GeneticMode mode)
                                    { return (bits & bit_for(mode)) != 0; });
    }

    // Derived from each() so undefined bits can never desync the two.
    [[nodiscard]] constexpr auto size() const -> std::size_t
    {
        return static_cast<std::size_t>(std::ranges::distance(each()));
    }

    // Inverse of each(): position of mode in the canonical traversal, or size()
    // if absent. Also derived from each(), so the two cannot disagree.
    [[nodiscard]] constexpr auto index_of(GeneticMode mode) const -> std::size_t
    {
        return static_cast<std::size_t>(std::ranges::distance(
            each()
            | std::views::take_while([mode](GeneticMode current)
                                     { return current != mode; })));
    }

    // Inverse of index_of: the mode occupying a position in the canonical
    // traversal. Also derived from each(), so the two cannot disagree.
    [[nodiscard]] constexpr auto at(std::size_t index) const -> GeneticMode
    {
        auto traversal = each();
        return *std::ranges::next(
            traversal.begin(), static_cast<std::ptrdiff_t>(index));
    }

    constexpr auto operator|=(GeneticModeSet other) noexcept -> GeneticModeSet&
    {
        bits = static_cast<std::uint8_t>(bits | other.bits);
        return *this;
    }

    [[nodiscard]] friend constexpr auto operator|(
        GeneticModeSet lhs,
        GeneticModeSet rhs) noexcept -> GeneticModeSet
    {
        lhs |= rhs;
        return lhs;
    }

    constexpr auto operator==(const GeneticModeSet&) const noexcept -> bool
        = default;

   private:
    static_assert(
        all_genetic_modes.size() <= std::numeric_limits<std::uint8_t>::digits);

    [[nodiscard]] static constexpr auto bit_for(GeneticMode mode) noexcept
        -> std::uint8_t
    {
        return static_cast<std::uint8_t>(1U << std::to_underlying(mode));
    }
};

[[nodiscard]] constexpr auto operator|(
    GeneticMode lhs,
    GeneticMode rhs) noexcept -> GeneticModeSet
{
    return GeneticModeSet{lhs} | GeneticModeSet{rhs};
}

// Canonical string tokens for a GeneticModeSet, including the composite "AD".
// Single source for both parsing (lexical_cast) and CLI validation.
inline constexpr std::array genetic_mode_set_names{
    std::pair{GeneticModeSet{GeneticMode::A}, std::string_view{"A"}},
    std::pair{GeneticModeSet{GeneticMode::D}, std::string_view{"D"}},
    std::pair{GeneticMode::A | GeneticMode::D, std::string_view{"AD"}},
};

// The canonical token of a set. Matching goes through each(), so undefined bits
// are ignored here as everywhere else.
[[nodiscard]] constexpr auto name_of(GeneticModeSet modes) -> std::string_view
{
    for (const auto& [value, name] : genetic_mode_set_names)
    {
        if (std::ranges::equal(value.each(), modes.each()))
        {
            return name;
        }
    }
    return {};
}

}  // namespace gelex

template <>
struct fmt::formatter<gelex::GeneticMode> : fmt::formatter<std::string_view>
{
    auto format(gelex::GeneticMode mode, auto& ctx) const
    {
        return fmt::formatter<std::string_view>::format(
            to_string_view(mode), ctx);
    }

   private:
    static constexpr auto to_string_view(gelex::GeneticMode mode)
        -> std::string_view
    {
        for (const auto& [value, name] : gelex::genetic_mode_names)
        {
            if (value == mode)
            {
                return name;
            }
        }
        return "unknown";
    }
};

template <>
struct fmt::formatter<gelex::GeneticModeSet> : fmt::formatter<std::string_view>
{
    auto format(gelex::GeneticModeSet modes, auto& ctx) const
    {
        return fmt::formatter<std::string_view>::format(
            gelex::name_of(modes), ctx);
    }
};

#endif  // GELEX_TYPES_GENETIC_MODE_H_
