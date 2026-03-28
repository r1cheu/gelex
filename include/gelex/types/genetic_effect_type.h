#ifndef GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
#define GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
#include <cstdint>
#include <format>
#include <string_view>

#include "gelex/exception.h"

namespace gelex
{

enum class GeneticKind : uint8_t
{
    NotGenetic = 255,
    Add = 0,
    Dom = 1
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

    Category category;
    GeneticKind genetic_kind{GeneticKind::NotGenetic};

    auto operator==(const EffectType&) const -> bool = default;

    static constexpr auto add() -> EffectType
    {
        return {Category::Genetic, GeneticKind::Add};
    }
    static constexpr auto dom() -> EffectType
    {
        return {Category::Genetic, GeneticKind::Dom};
    }
    static constexpr auto fixed() -> EffectType { return {Category::Fixed}; }
    static constexpr auto random() -> EffectType { return {Category::Random}; }
    static constexpr auto residual() -> EffectType
    {
        return {Category::Residual};
    }
    static constexpr auto from_genetic(GeneticKind gk) -> EffectType
    {
        return {Category::Genetic, gk};
    }

    constexpr auto to_byte() const -> uint8_t
    {
        if (category == Category::Genetic)
        {
            return static_cast<uint8_t>(genetic_kind);
        }
        return static_cast<uint8_t>(category);
    }

    static constexpr auto from_byte(uint8_t wire) -> EffectType
    {
        switch (wire)
        {
            case 0:
                return add();
            case 1:
                return dom();
            case 2:
                return fixed();
            case 3:
                return random();
            case 4:
                return residual();
            default:
                throw GelexException(
                    "EffectType::from_byte: unknown wire value");
        }
    }
};

namespace genetic_kind
{

inline auto to_string(GeneticKind type) -> std::string_view
{
    switch (type)
    {
        case GeneticKind::Add:
            return "Additive";
        case GeneticKind::Dom:
            return "Dominance";
        default:
            return "Unknown";
    }
}

inline auto to_file_suffix(GeneticKind type) -> std::string_view
{
    switch (type)
    {
        case GeneticKind::Add:
            return "add";
        case GeneticKind::Dom:
            return "dom";
        default:
            return "unk";
    }
}

inline auto to_variance_label(GeneticKind type) -> std::string_view
{
    switch (type)
    {
        case GeneticKind::Add:
            return "σ²_add";
        case GeneticKind::Dom:
            return "σ²_dom";
        default:
            return "σ²_unk";
    }
}

inline auto to_heritability_label(GeneticKind type) -> std::string_view
{
    switch (type)
    {
        case GeneticKind::Add:
            return "h²";
        case GeneticKind::Dom:
            return "δ²";
        default:
            return "?²";
    }
}

}  // namespace genetic_kind

}  // namespace gelex

template <>
struct std::formatter<gelex::EffectType> : std::formatter<std::string_view>
{
    auto format(const gelex::EffectType& et, auto& ctx) const
    {
        using Cat = gelex::EffectType::Category;
        std::string_view name;
        switch (et.category)
        {
            case Cat::Genetic:
                name = (et.genetic_kind == gelex::GeneticKind::Add)
                           ? "Additive"
                           : "Dominance";
                break;
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
        return std::formatter<std::string_view>::format(name, ctx);
    }
};

#endif  // GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
