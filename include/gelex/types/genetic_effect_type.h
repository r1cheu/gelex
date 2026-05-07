#ifndef GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
#define GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
#include <algorithm>
#include <array>
#include <cstdint>
#include <optional>
#include <string_view>
#include <utility>

#include <fmt/format.h>

namespace gelex
{

enum class GeneticMode : uint8_t
{
    A,
    D,
    AD,
    kCount,
};

inline constexpr std::array<
    std::pair<GeneticMode, std::string_view>,
    std::to_underlying(GeneticMode::kCount)>
    kGeneticModeNames = {{
        {GeneticMode::A, "A"},
        {GeneticMode::D, "D"},
        {GeneticMode::AD, "AD"},
    }};

inline auto get_genetic_mode(std::string_view sv) -> std::optional<GeneticMode>
{
    const auto* it = std::ranges::find_if(
        kGeneticModeNames, [sv](const auto& p) { return p.second == sv; });
    if (it != kGeneticModeNames.end())
    {
        return it->first;
    }
    return std::nullopt;
}

struct EffectType
{
    enum class Category : uint8_t
    {
        Genetic = 0,
        Fixed = 2,
        Random = 3,
        Residual = 4
    };

    Category category{};
    std::optional<GeneticMode> genetic_mode;

    auto operator==(const EffectType&) const -> bool = default;

    static constexpr auto add() -> EffectType
    {
        return {Category::Genetic, GeneticMode::A};
    }
    static constexpr auto dom() -> EffectType
    {
        return {Category::Genetic, GeneticMode::D};
    }
    static constexpr auto fixed() -> EffectType
    {
        return {Category::Fixed, {}};
    }
    static constexpr auto random() -> EffectType
    {
        return {Category::Random, {}};
    }
    static constexpr auto residual() -> EffectType
    {
        return {Category::Residual, {}};
    }
    static constexpr auto from_genetic(GeneticMode gm) -> EffectType
    {
        return {Category::Genetic, gm};
    }
};

inline constexpr std::array kAllGeneticModes{GeneticMode::A, GeneticMode::D};

namespace genetic_mode
{

inline auto to_file_suffix(GeneticMode type) -> std::string_view
{
    switch (type)
    {
        case GeneticMode::A:
            return "add";
        case GeneticMode::D:
            return "dom";
        case GeneticMode::AD:
            return "ad";
        case GeneticMode::kCount:
            break;
    }
    return "unknown";
}

inline auto to_variance_label(GeneticMode type) -> std::string_view
{
    switch (type)
    {
        case GeneticMode::A:
            return "σ²_add";
        case GeneticMode::D:
            return "σ²_dom";
        case GeneticMode::AD:
            return "σ²_g";
        case GeneticMode::kCount:
            break;
    }
    return "unknown";
}

inline auto to_heritability_label(GeneticMode type) -> std::string_view
{
    switch (type)
    {
        case GeneticMode::A:
            return "h²";
        case GeneticMode::D:
            return "δ²";
        case GeneticMode::AD:
            return "H²";
        case GeneticMode::kCount:
            break;
    }
    return "unknown";
}

}  // namespace genetic_mode

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
        for (const auto& [value, name] : gelex::kGeneticModeNames)
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
struct fmt::formatter<gelex::EffectType> : fmt::formatter<std::string_view>
{
    auto format(const gelex::EffectType& et, auto& ctx) const
    {
        using Cat = gelex::EffectType::Category;
        if (et.category == Cat::Genetic)
        {
            return fmt::format_to(ctx.out(), "{}", *et.genetic_mode);
        }
        return fmt::formatter<std::string_view>::format(
            to_string_view(et.category), ctx);
    }

   private:
    static constexpr auto to_string_view(gelex::EffectType::Category cat)
        -> std::string_view
    {
        using Cat = gelex::EffectType::Category;
        switch (cat)
        {
            case Cat::Genetic:
                return "Genetic";
            case Cat::Fixed:
                return "Fixed";
            case Cat::Random:
                return "Random";
            case Cat::Residual:
                return "Residual";
        }
        return "unknown";
    }
};

#endif  // GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
