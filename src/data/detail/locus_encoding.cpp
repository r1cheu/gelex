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

#include "gelex/data/locus_encoding.h"

#include <array>
#include <cmath>

namespace gelex
{

auto encoding_spec_from_method(GeneticMode effect, GenotypeMethod method)
    -> EncodingSpec
{
    EncodingSpec spec;
    spec.effect = effect;
    spec.normalization = is_center(method) ? Normalization::Center
                                           : Normalization::CenterScale;
    spec.moment_basis
        = is_hwe(method) ? MomentBasis::Theoretical : MomentBasis::Empirical;

    if (is_noia(method))
    {
        spec.dominance_code = DominanceCode::NOIA;
    }
    else if (is_orthogonal(method))
    {
        spec.dominance_code = DominanceCode::HWE;
    }
    else
    {
        spec.dominance_code = DominanceCode::Het;
    }

    return spec;
}

auto LocusStats::pA2A2() const -> double
{
    return static_cast<double>(nA2A2) / static_cast<double>(n_nonmissing());
}

auto LocusStats::pA1A2() const -> double
{
    return static_cast<double>(nA1A2) / static_cast<double>(n_nonmissing());
}

auto LocusStats::pA1A1() const -> double
{
    return static_cast<double>(nA1A1) / static_cast<double>(n_nonmissing());
}

auto LocusStats::A1freq() const -> double
{
    return pA1A1() + (0.5 * pA1A2());
}

namespace detail
{

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

auto weighted_mean(
    const std::array<double, 3>& values,
    const MomentWeights& weights) -> double
{
    return (weights.A2A2 * values[0]) + (weights.A1A2 * values[1])
           + (weights.A1A1 * values[2]);
}

auto weighted_var(
    const std::array<double, 3>& values,
    const MomentWeights& weights,
    double mean) -> double
{
    const double x0{values[0] - mean};
    const double x1{values[1] - mean};
    const double x2{values[2] - mean};

    return (weights.A2A2 * x0 * x0) + (weights.A1A2 * x1 * x1)
           + (weights.A1A1 * x2 * x2);
}

auto make_dominance_het() -> CodeMap
{
    return CodeMap{.value = {0.0, 1.0, 0.0}, .valid = true};
}

auto make_dominance_hwe(const LocusStats& stats) -> CodeMap
{
    const double p{stats.A1freq()};

    return CodeMap{.value = {0.0, 2.0 * p, (4.0 * p) - 2.0}, .valid = true};
}

auto make_dominance_noia(const LocusStats& stats, double tol) -> CodeMap
{
    CodeMap out;

    const double pA2A2{stats.pA2A2()};
    const double pA1A2{stats.pA1A2()};
    const double pA1A1{stats.pA1A1()};

    const double diff{pA1A1 - pA2A2};
    const double denom{pA1A1 + pA2A2 - (diff * diff)};

    if (denom < tol)
    {
        out.valid = false;
        out.value = {0.0, 0.0, 0.0};
        return out;
    }

    const double cA1A1{-2.0 * pA2A2 * pA1A2 / denom};
    const double cA1A2{4.0 * pA1A1 * pA2A2 / denom};
    const double cA2A2{-2.0 * pA1A1 * pA1A2 / denom};

    out.value = {cA2A2, cA1A2, cA1A1};
    return out;
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

    CodeMap code;
    if (spec.effect == gelex::GeneticMode::A)
    {
        code.value = {0.0, 1.0, 2.0};
    }
    else
    {
        switch (spec.dominance_code)
        {
            case DominanceCode::Het:
                code = make_dominance_het();
                break;
            case DominanceCode::HWE:
                code = make_dominance_hwe(stats);
                break;
            case DominanceCode::NOIA:
                code = make_dominance_noia(stats, tol);
                break;
        }
    }

    if (!code.valid)
    {
        out.valid = false;
        out.sd = 0.0;
        return out;
    }

    const MomentWeights weights{make_moment_weights(stats, spec.moment_basis)};

    out.mean = weighted_mean(code.value, weights);
    out.var = weighted_var(code.value, weights, out.mean);

    if (out.var < tol)
    {
        out.valid = false;
        out.code = {0.0, 0.0, 0.0};
        out.sd = 0.0;
        return out;
    }

    out.sd = std::sqrt(out.var);
    out.code = code.value;

    if (spec.normalization == Normalization::Center
        || spec.normalization == Normalization::CenterScale)
    {
        for (double& value : out.code)
        {
            value -= out.mean;
        }
    }

    if (spec.normalization == Normalization::CenterScale)
    {
        for (double& value : out.code)
        {
            value /= out.sd;
        }
    }

    if (spec.normalization == Normalization::None)
    {
        out.missing_encoded_value = out.mean;
    }
    else
    {
        out.missing_encoded_value = 0.0;
    }

    out.valid = true;
    return out;
}

}  // namespace detail

}  // namespace gelex
