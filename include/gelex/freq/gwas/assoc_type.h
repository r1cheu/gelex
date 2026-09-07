// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_GWAS_ASSOC_TYPE_H_
#define GELEX_FREQ_GWAS_ASSOC_TYPE_H_

#include <cstdint>

namespace gelex
{

enum class AssocType : uint8_t
{
    Single = 0,
    Joint = 1
};

}  // namespace gelex

#endif  // GELEX_FREQ_GWAS_ASSOC_TYPE_H_
