// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_ENCODE_TYPES_H_
#define GELEX_DATA_ENCODE_TYPES_H_

#include <Eigen/Core>
#include <cstdint>

#include "gelex/data/encode/stats.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/snp_lut.h"
#include "gelex/genetic_mode.h"

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

    SnpLut lut{0.0, 0.0, 0.0, 0.0};

    double mean{0};
    double var{0};
    double sd{1};

    bool valid{false};
};

}  // namespace gelex

#endif  // GELEX_DATA_ENCODE_TYPES_H_
