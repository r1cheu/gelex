#ifndef GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
#define GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
#include <cstdint>
#include <string_view>

namespace gelex
{
enum class GeneticEffectType : uint8_t
{
    Add,
    Dom
};

namespace genetic_effect_type
{

inline auto to_string(GeneticEffectType type) -> std::string_view
{
    switch (type)
    {
        case GeneticEffectType::Add:
            return "Additive";
        case GeneticEffectType::Dom:
            return "Dominance";
    }
    return "Unknown";
}

inline auto to_file_suffix(GeneticEffectType type) -> std::string_view
{
    switch (type)
    {
        case GeneticEffectType::Add:
            return "add";
        case GeneticEffectType::Dom:
            return "dom";
    }
    return "unk";
}

inline auto to_variance_label(GeneticEffectType type) -> std::string_view
{
    switch (type)
    {
        case GeneticEffectType::Add:
            return "σ²_add";
        case GeneticEffectType::Dom:
            return "σ²_dom";
    }
    return "σ²_unk";
}

inline auto to_heritability_label(GeneticEffectType type) -> std::string_view
{
    switch (type)
    {
        case GeneticEffectType::Add:
            return "h²";
        case GeneticEffectType::Dom:
            return "δ²";
    }
    return "?²";
}

}  // namespace genetic_effect_type

enum class ModelType : uint8_t
{
    A,
    D,
    AD
};

struct LocusStatistic
{
    double mean{0};
    double stddev{0};
    double maf{0};
    bool is_monomorphic{false};
};

}  // namespace gelex

#endif  // GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
