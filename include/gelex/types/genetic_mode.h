#ifndef GELEX_TYPES_GENETIC_MODE_H_
#define GELEX_TYPES_GENETIC_MODE_H_
#include <array>
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <map>
#include <ranges>
#include <string_view>
#include <utility>

namespace gelex
{

enum class GeneticMode : uint8_t
{
    A,
    D,
};

inline constexpr std::array GENETIC_MODE_NAMES{
    std::pair{GeneticMode::A, std::string_view{"A"}},
    std::pair{GeneticMode::D, std::string_view{"D"}},
};

inline constexpr std::array ALL_GENETIC_MODES{GeneticMode::A, GeneticMode::D};

template <typename T>
using ModeMap = std::map<GeneticMode, T>;

// A subset of GeneticMode, backed by std::bitset so {A}, {D}, {A, D} are the
// only meaningful values. each() yields its members in enum order.
class GeneticModeSet
{
   public:
    GeneticModeSet() = default;

    explicit constexpr GeneticModeSet(GeneticMode mode)
    {
        bits_.set(std::to_underlying(mode));
    }

    [[nodiscard]] auto contains(GeneticMode mode) const -> bool
    {
        return bits_.test(std::to_underlying(mode));
    }

    [[nodiscard]] auto size() const -> std::size_t { return bits_.count(); }

    constexpr auto operator|=(GeneticModeSet other) -> GeneticModeSet&
    {
        bits_ |= other.bits_;
        return *this;
    }

    [[nodiscard]] friend constexpr auto operator|(
        GeneticModeSet lhs,
        GeneticModeSet rhs) -> GeneticModeSet
    {
        lhs |= rhs;
        return lhs;
    }

    auto operator==(const GeneticModeSet&) const -> bool = default;

    [[nodiscard]] auto each() const
    {
        return ALL_GENETIC_MODES
               | std::views::filter(
                   [bits = bits_](GeneticMode mode)
                   { return bits.test(std::to_underlying(mode)); });
    }

   private:
    std::bitset<ALL_GENETIC_MODES.size()> bits_;
};

[[nodiscard]] constexpr auto operator|(GeneticMode lhs, GeneticMode rhs)
    -> GeneticModeSet
{
    return GeneticModeSet{lhs} | GeneticModeSet{rhs};
}

// Canonical string tokens for a GeneticModeSet, including the composite "AD".
// Single source for both parsing (lexical_cast) and CLI validation.
inline constexpr std::array GENETIC_MODE_SET_NAMES{
    std::pair{GeneticModeSet{GeneticMode::A}, std::string_view{"A"}},
    std::pair{GeneticModeSet{GeneticMode::D}, std::string_view{"D"}},
    std::pair{GeneticMode::A | GeneticMode::D, std::string_view{"AD"}},
};

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
        for (const auto& [value, name] : gelex::GENETIC_MODE_NAMES)
        {
            if (value == mode)
            {
                return name;
            }
        }
        return "unknown";
    }
};

#endif  // GELEX_TYPES_GENETIC_MODE_H_
