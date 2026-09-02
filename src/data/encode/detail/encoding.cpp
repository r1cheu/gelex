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

#include "gelex/data/encode/detail/encoding.h"

#include <array>
#include <cmath>
#include <optional>

#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

namespace
{

constexpr Eigen::Index missing_code{1};
constexpr std::array<Eigen::Index, 3> non_missing_codes{0, 2, 3};

}  // namespace

auto make_moment_weights(const LocusStats& stats, MomentBasis basis)
    -> MomentWeights
{
    if (basis == MomentBasis::Empirical)
    {
        return {
            .A2A2 = stats.pA2A2(),
            .A1A2 = stats.pA1A2(),
            .A1A1 = stats.pA1A1()};
    }

    const double p{stats.A1freq()};
    const double q{1.0 - p};

    return {.A2A2 = q * q, .A1A2 = 2.0 * p * q, .A1A1 = p * p};
}

auto weighted_mean(const SnpLut& lut, const MomentWeights& weights) -> double
{
    return (weights.A1A1 * lut[0]) + (weights.A1A2 * lut[2])
           + (weights.A2A2 * lut[3]);
}

auto weighted_var(const SnpLut& lut, const MomentWeights& weights, double mean)
    -> double
{
    const double x0{lut[0] - mean};
    const double x2{lut[2] - mean};
    const double x3{lut[3] - mean};

    return (weights.A1A1 * x0 * x0) + (weights.A1A2 * x2 * x2)
           + (weights.A2A2 * x3 * x3);
}

auto make_dominance_het() -> SnpLut
{
    return {0.0, 0.0, 1.0, 0.0};
}

auto make_dominance_hwe(const LocusStats& stats) -> SnpLut
{
    const double p{stats.A1freq()};

    return {(4.0 * p) - 2.0, 0.0, 2.0 * p, 0.0};
}

auto make_dominance_noia(const LocusStats& stats, double tol)
    -> std::optional<SnpLut>
{
    const double pA2A2{stats.pA2A2()};
    const double pA1A2{stats.pA1A2()};
    const double pA1A1{stats.pA1A1()};

    const double diff{pA1A1 - pA2A2};
    const double denom{pA1A1 + pA2A2 - (diff * diff)};

    if (denom < tol)
    {
        return std::nullopt;
    }

    const double cA1A1{-2.0 * pA2A2 * pA1A2 / denom};
    const double cA1A2{4.0 * pA1A1 * pA2A2 / denom};
    const double cA2A2{-2.0 * pA1A1 * pA1A2 / denom};

    return SnpLut{cA1A1, 0.0, cA1A2, cA2A2};
}

auto make_locus_encoding(
    Eigen::Index marker_index,
    const LocusStats& stats,
    const EncodingSpec& spec,
    double tol) -> LocusEncoding
{
    LocusEncoding out;
    out.column_index = marker_index;
    out.marker_index = marker_index;
    out.stats = stats;

    if (!stats.has_nonmissing())
    {
        out.valid = false;
        out.sd = 0.0;
        return out;
    }

    std::optional<SnpLut> lut;
    if (spec.effect == gelex::GeneticMode::A)
    {
        lut = SnpLut{2.0, 0.0, 1.0, 0.0};
    }
    else
    {
        switch (spec.dominance_code)
        {
            case DominanceCode::Het:
                lut = make_dominance_het();
                break;
            case DominanceCode::HWE:
                lut = make_dominance_hwe(stats);
                break;
            case DominanceCode::NOIA:
                lut = make_dominance_noia(stats, tol);
                break;
        }
    }

    if (!lut)
    {
        out.valid = false;
        out.sd = 0.0;
        return out;
    }

    const MomentWeights weights{make_moment_weights(stats, spec.moment_basis)};

    out.mean = weighted_mean(*lut, weights);
    out.var = weighted_var(*lut, weights, out.mean);

    if (out.var < tol)
    {
        out.valid = false;
        out.sd = 0.0;
        return out;
    }

    out.sd = std::sqrt(out.var);
    out.lut = *lut;

    switch (spec.normalization)
    {
        case Normalization::Center:
            for (const Eigen::Index raw_code : non_missing_codes)
            {
                out.lut[raw_code] -= out.mean;
            }
            out.lut[missing_code] = 0.0;
            break;
        case Normalization::CenterScale:
            for (const Eigen::Index raw_code : non_missing_codes)
            {
                out.lut[raw_code] = (out.lut[raw_code] - out.mean) / out.sd;
            }
            out.lut[missing_code] = 0.0;
            break;
        case Normalization::None:
            out.lut[missing_code] = out.mean;
            break;
        default:
            out.lut[missing_code] = 0.0;
            break;
    }
    out.valid = true;
    return out;
}

}  // namespace gelex::detail
