// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_DRAW_TRAITS_H_
#define GELEX_BAYES_GENETIC_DRAW_TRAITS_H_

#include <type_traits>
#include <variant>

#include "gelex/bayes/genetic/types.h"
#include "gelex/io/binary_format.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

template <detail::SupportedDtype>
class DenseStream;

template <VarianceLayout Kind>
using marker_variance_dtype_t
    = std::conditional_t<Kind == VarianceLayout::Pooled, double, float>;
template <VarianceLayout Kind>
using marker_variance_writer_t = DenseStream<marker_variance_dtype_t<Kind>>;

template <MixtureWeightUpdate Update>
using probability_writer_t = std::conditional_t<
    Update == MixtureWeightUpdate::Enabled,
    DenseStream<double>,
    std::monostate>;

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_DRAW_TRAITS_H_
