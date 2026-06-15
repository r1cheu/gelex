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

#ifndef GELEX_DATA_DETAIL_LOCUS_ENCODING_H_
#define GELEX_DATA_DETAIL_LOCUS_ENCODING_H_

#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>

#include <Eigen/Core>

#include "gelex/data/locus_encoding_types.h"

namespace gelex::detail
{

struct MomentWeights
{
    double A2A2{0};
    double A1A2{0};
    double A1A1{0};
};

struct CodeMap
{
    std::array<double, 3> value{0.0, 0.0, 0.0};
    bool valid{true};
};

[[nodiscard]] auto make_moment_weights(
    const LocusStats& stats,
    MomentBasis basis) -> MomentWeights;

[[nodiscard]] auto weighted_mean(
    const std::array<double, 3>& values,
    const MomentWeights& weights) -> double;

[[nodiscard]] auto weighted_var(
    const std::array<double, 3>& values,
    const MomentWeights& weights,
    double mean) -> double;

[[nodiscard]] auto make_dominance_het() -> CodeMap;
[[nodiscard]] auto make_dominance_hwe(const LocusStats& stats) -> CodeMap;
[[nodiscard]] auto make_dominance_noia(
    const LocusStats& stats,
    double tol = 1e-12) -> CodeMap;

[[nodiscard]] auto make_locus_encoding(
    Eigen::Index marker_index,
    const LocusStats& stats,
    const EncodingSpec& spec,
    double tol = 1e-12) -> LocusEncoding;

template <std::floating_point T>
auto make_loci_encoding(
    const Eigen::Ref<const Eigen::MatrixX<T>>& genotypes,
    const EncodingSpec& spec,
    T tol = static_cast<T>(1e-12),
    Eigen::Index marker_offset = 0) -> LociEncoding
{
    LociEncoding out;
    out.spec = spec;

    const Eigen::Index num_markers{genotypes.cols()};
    out.loci.reserve(static_cast<std::size_t>(num_markers));

    for (Eigen::Index marker_index{0}; marker_index < num_markers;
         ++marker_index)
    {
        const auto locus = genotypes.col(marker_index);
        const LocusStats stats{compute_locus_stats<T>(locus)};
        const LocusEncoding encoding{make_locus_encoding(
            marker_offset + marker_index,
            stats,
            spec,
            static_cast<double>(tol))};

        out.loci.push_back(encoding);
        out.loci.back().column_index = marker_index;
    }

    return out;
}

template <std::floating_point InputT, std::floating_point OutputT>
auto transform_into(
    const Eigen::Ref<const Eigen::MatrixX<InputT>>& genotypes,
    Eigen::Ref<Eigen::MatrixX<OutputT>> output,
    const LociEncoding& encoding) -> void
{
    assert(genotypes.rows() == output.rows());
    assert(genotypes.cols() == output.cols());

    for (const auto& locus_encoding : encoding.loci)
    {
        const Eigen::Index column_index{locus_encoding.column_index};
        assert(column_index >= 0);
        assert(column_index < genotypes.cols());

        if (!locus_encoding.valid)
        {
            output.col(column_index).setZero();
            continue;
        }

        for (Eigen::Index sample_index{0}; sample_index < genotypes.rows();
             ++sample_index)
        {
            const InputT genotype{genotypes(sample_index, column_index)};
            auto encoded_value{
                static_cast<OutputT>(locus_encoding.missing_encoded_value)};

            if (genotype == InputT{0})
            {
                encoded_value = static_cast<OutputT>(locus_encoding.code[0]);
            }
            else if (genotype == InputT{1})
            {
                encoded_value = static_cast<OutputT>(locus_encoding.code[1]);
            }
            else if (genotype == InputT{2})
            {
                encoded_value = static_cast<OutputT>(locus_encoding.code[2]);
            }

            output(sample_index, column_index) = encoded_value;
        }
    }
}

template <std::floating_point T>
auto transform_inplace(
    Eigen::Ref<Eigen::MatrixX<T>> genotypes,
    const LociEncoding& encoding) -> void
{
    for (const auto& locus_encoding : encoding.loci)
    {
        const Eigen::Index column_index{locus_encoding.column_index};
        assert(column_index >= 0);
        assert(column_index < genotypes.cols());

        if (!locus_encoding.valid)
        {
            genotypes.col(column_index).setZero();
            continue;
        }

        for (Eigen::Index sample_index{0}; sample_index < genotypes.rows();
             ++sample_index)
        {
            const T genotype{genotypes(sample_index, column_index)};
            auto encoded_value{
                static_cast<T>(locus_encoding.missing_encoded_value)};

            if (genotype == T{0})
            {
                encoded_value = static_cast<T>(locus_encoding.code[0]);
            }
            else if (genotype == T{1})
            {
                encoded_value = static_cast<T>(locus_encoding.code[1]);
            }
            else if (genotype == T{2})
            {
                encoded_value = static_cast<T>(locus_encoding.code[2]);
            }

            genotypes(sample_index, column_index) = encoded_value;
        }
    }
}

}  // namespace gelex::detail

#endif  // GELEX_DATA_DETAIL_LOCUS_ENCODING_H_
