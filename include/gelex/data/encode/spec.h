// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_ENCODE_SPEC_H_
#define GELEX_DATA_ENCODE_SPEC_H_

#include "gelex/data/encode/types.h"
#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

[[nodiscard]] auto encoding_spec_from_method(
    GeneticMode effect,
    GenotypeMethod method) -> EncodingSpec;

}  // namespace gelex

#endif  // GELEX_DATA_ENCODE_SPEC_H_
