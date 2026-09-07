// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_ENCODE_DETAIL_ENCODING_H_
#define GELEX_DATA_ENCODE_DETAIL_ENCODING_H_

#include <Eigen/Core>
#include <optional>

#include "gelex/data/encode/types.h"
#include "gelex/data/snp_lut.h"

namespace gelex::detail
{

struct MomentWeights
{
    double A2A2{0};
    double A1A2{0};
    double A1A1{0};
};

[[nodiscard]] auto make_moment_weights(
    const LocusStats& stats,
    MomentBasis basis) -> MomentWeights;

[[nodiscard]] auto weighted_mean(
    const SnpLut& lut,
    const MomentWeights& weights) -> double;

[[nodiscard]] auto weighted_var(
    const SnpLut& lut,
    const MomentWeights& weights,
    double mean) -> double;

[[nodiscard]] auto make_dominance_het() -> SnpLut;
[[nodiscard]] auto make_dominance_hwe(const LocusStats& stats) -> SnpLut;
[[nodiscard]] auto make_dominance_noia(
    const LocusStats& stats,
    double tol = 1e-12) -> std::optional<SnpLut>;

[[nodiscard]] auto make_locus_encoding(
    Eigen::Index marker_index,
    const LocusStats& stats,
    const EncodingSpec& spec,
    double tol = 1e-12) -> LocusEncoding;

}  // namespace gelex::detail

#endif  // GELEX_DATA_ENCODE_DETAIL_ENCODING_H_
