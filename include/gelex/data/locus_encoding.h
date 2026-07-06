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

#ifndef GELEX_DATA_LOCUS_ENCODING_H_
#define GELEX_DATA_LOCUS_ENCODING_H_

#include <cassert>
#include <concepts>

#include <Eigen/Core>

#include "gelex/data/genotype_method.h"
#include "gelex/data/locus_encoding_types.h"
#include "gelex/types/genetic_mode.h"

#include "gelex/data/detail/locus_encoding.h"

namespace gelex
{

[[nodiscard]] auto encoding_spec_from_method(
    GeneticMode effect,
    GenotypeMethod method) -> EncodingSpec;

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

template <std::floating_point T>
auto encode_inplace(
    Eigen::Ref<Eigen::MatrixX<T>> genotypes,
    GeneticMode effect,
    GenotypeMethod method,
    T tol = static_cast<T>(1e-12),
    Eigen::Index marker_offset = 0) -> LociEncoding
{
    const EncodingSpec spec{encoding_spec_from_method(effect, method)};
    LociEncoding encoding{
        detail::make_loci_encoding<T>(genotypes, spec, tol, marker_offset)};
    transform_inplace<T>(genotypes, encoding);
    return encoding;
}

template <std::floating_point InputT, std::floating_point OutputT>
auto encode_into(
    const Eigen::Ref<const Eigen::MatrixX<InputT>>& genotypes,
    Eigen::Ref<Eigen::MatrixX<OutputT>> output,
    GeneticMode effect,
    GenotypeMethod method,
    InputT tol = static_cast<InputT>(1e-12),
    Eigen::Index marker_offset = 0) -> LociEncoding
{
    const EncodingSpec spec{encoding_spec_from_method(effect, method)};
    LociEncoding encoding{detail::make_loci_encoding<InputT>(
        genotypes, spec, tol, marker_offset)};
    transform_into<InputT, OutputT>(genotypes, output, encoding);
    return encoding;
}

}  // namespace gelex

#endif  // GELEX_DATA_LOCUS_ENCODING_H_
