#ifndef GELEX_TYPES_GENETIC_MODE_H_
#define GELEX_TYPES_GENETIC_MODE_H_
#include <array>
#include <cstdint>
#include <string_view>

#include <fmt/format.h>

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
