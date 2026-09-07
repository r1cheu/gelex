// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_TYPES_H_
#define GELEX_BAYES_GENETIC_TYPES_H_

#include <cstddef>
#include <cstdint>

namespace gelex
{

enum class VarianceLayout : std::uint8_t
{
    Pooled,
    Unpooled,
};

enum class MixtureWeightUpdate : std::uint8_t
{
    Disabled,
    Enabled,
};

struct GeneticDimensions
{
    std::size_t individual;
    std::size_t marker;
};

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_TYPES_H_
