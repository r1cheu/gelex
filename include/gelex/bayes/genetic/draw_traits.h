// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_DRAW_TRAITS_H_
#define GELEX_BAYES_GENETIC_DRAW_TRAITS_H_

#include <cstdint>
#include <string_view>
#include <type_traits>
#include <variant>

#include "gelex/bayes/genetic/types.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/csc_writer.h"
#include "gelex/io/dense_writer.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

// Draw factories reserve dense payloads in `dense` and marker-level sparse
// payloads in `sparse`; the resulting streams borrow both writers.
struct DrawWriters
{
    DenseWriter& dense;
    CscWriter& sparse;
};

template <VarianceLayout Kind>
using marker_variance_dtype_t
    = std::conditional_t<Kind == VarianceLayout::Pooled, double, float>;
template <VarianceLayout Kind>
using marker_variance_writer_t = DenseStream<marker_variance_dtype_t<Kind>>;

template <CoefficientLayout Layout>
using coefficients_writer_t = std::conditional_t<
    Layout == CoefficientLayout::Dense,
    DenseStream<float>,
    CscStream<float>>;

using assignment_writer_t = CscStream<std::uint8_t>;

template <MixtureWeightUpdate Update>
using probability_writer_t = std::conditional_t<
    Update == MixtureWeightUpdate::Enabled,
    DenseStream<double>,
    std::monostate>;

template <CoefficientLayout Layout>
auto reserve_coefficients(
    DrawWriters writers,
    std::string_view identifier,
    BinaryShape shape) -> coefficients_writer_t<Layout>
{
    if constexpr (Layout == CoefficientLayout::Dense)
    {
        return writers.dense.reserve<float>(identifier, shape);
    }
    else
    {
        return writers.sparse.reserve<float>(identifier, shape);
    }
}

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_DRAW_TRAITS_H_
