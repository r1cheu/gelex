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
