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

// Only the pooled marker variance is a model parameter worth keeping;
// unpooled variances are marker-level nuisance draws and are not stored.
template <VarianceLayout Kind>
using marker_variance_writer_t = std::conditional_t<
    Kind == VarianceLayout::Pooled,
    DenseStream<double>,
    std::monostate>;

template <VarianceLayout Kind>
auto reserve_marker_variance(
    DrawWriters writers,
    std::string_view identifier,
    std::uint64_t draw_count) -> marker_variance_writer_t<Kind>
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return writers.dense.reserve<double>(
            identifier, BinaryShape{1, draw_count});
    }
    else
    {
        return {};
    }
}

template <CoefficientLayout Layout>
using coefficients_writer_t = std::conditional_t<
    Layout == CoefficientLayout::Dense,
    DenseStream<double>,
    CscStream<double>>;

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
        return writers.dense.reserve<double>(identifier, shape);
    }
    else
    {
        return writers.sparse.reserve<double>(identifier, shape);
    }
}

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_DRAW_TRAITS_H_
