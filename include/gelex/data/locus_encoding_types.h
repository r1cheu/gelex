/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GELEX_DATA_LOCUS_ENCODING_TYPES_H_
#define GELEX_DATA_LOCUS_ENCODING_TYPES_H_

#include <array>
#include <cstdint>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/genotype_method.h"
#include "gelex/data/locus_stats.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

enum class DominanceCode : std::uint8_t
{
    Het,
    HWE,
    NOIA
};

enum class Normalization : std::uint8_t
{
    None,
    Center,
    CenterScale
};

enum class MomentBasis : std::uint8_t
{
    Empirical,
    Theoretical
};

struct EncodingSpec
{
    GeneticMode effect{GeneticMode::A};

    DominanceCode dominance_code{DominanceCode::NOIA};

    Normalization normalization{Normalization::CenterScale};
    MomentBasis moment_basis{MomentBasis::Empirical};
};

struct LocusEncoding
{
    Eigen::Index column_index{-1};
    Eigen::Index marker_index{-1};

    LocusStats stats;

    std::array<double, 3> code{0.0, 0.0, 0.0};

    double mean{0};
    double var{0};
    double sd{1};

    double missing_encoded_value{0};

    bool valid{false};
};

struct LociEncoding
{
    EncodingSpec spec;
    std::vector<LocusEncoding> loci;
};

}  // namespace gelex

#endif  // GELEX_DATA_LOCUS_ENCODING_TYPES_H_
