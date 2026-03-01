#ifndef GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
#define GELEX_TYPES_GENETIC_EFFECT_TYPE_H_
#include <cstdint>

namespace gelex
{
enum class GeneticEffectType : uint8_t
{
    Add,
    Dom
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
