#ifndef GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
#define GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
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
    A = 0,
    D = 1,
    AD = 2
};

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
    }
    std::unreachable();
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
    }
    std::unreachable();
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
    }
    std::unreachable();
}

}  // namespace genetic_mode

}  // namespace gelex

template <>
struct fmt::formatter<gelex::GeneticMode> : fmt::formatter<std::string_view>
{
    auto format(gelex::GeneticMode mode, auto& ctx) const
    {
        std::string_view name;
        switch (mode)
        {
            case gelex::GeneticMode::A:
                name = "Additive";
                break;
            case gelex::GeneticMode::D:
                name = "Dominance";
                break;
            case gelex::GeneticMode::AD:
                name = "Additive Dominance";
                break;
        }
        return fmt::formatter<std::string_view>::format(name, ctx);
    }
};

template <>
struct fmt::formatter<gelex::EffectType> : fmt::formatter<std::string_view>
{
    auto format(const gelex::EffectType& et, auto& ctx) const
    {
        using Cat = gelex::EffectType::Category;
        std::string_view name;
        switch (et.category)
        {
            case Cat::Genetic:
                return fmt::format_to(ctx.out(), "{}", *et.genetic_mode);
            case Cat::Fixed:
                name = "Fixed";
                break;
            case Cat::Random:
                name = "Random";
                break;
            case Cat::Residual:
                name = "Residual";
                break;
        }
        return fmt::formatter<std::string_view>::format(name, ctx);
    }
};

#endif  // GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
