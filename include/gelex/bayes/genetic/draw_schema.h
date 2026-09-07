// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_DRAW_SCHEMA_H_
#define GELEX_BAYES_GENETIC_DRAW_SCHEMA_H_

/* Shared serialization schema for genetic draws: payload IDs,
 * storage dtypes, and corresponding writer types.
 */

#include <fmt/compile.h>
#include <string_view>
#include <type_traits>
#include <variant>

#include "gelex/bayes/genetic/types.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

template <detail::SupportedDtype>
class PayloadWriter;

GELEX_NAMESPACE_BEGIN(detail)
template <GeneticMode Mode>
inline constexpr auto genetic_id_storage
    = FMT_STATIC_FORMAT("genetic/{}", Mode);
GELEX_NAMESPACE_END(detail)

template <GeneticMode Mode>
inline const std::string_view genetic_id{
    detail::genetic_id_storage<Mode>.c_str()};

inline constexpr std::string_view joint_genetic_id = "genetic/joint";
inline constexpr std::string_view coefficients_id = "coefficients";
inline constexpr std::string_view variance_id = "variance";
inline constexpr std::string_view assignment_id = "assignment";
inline constexpr std::string_view probability_id = "probability";
inline constexpr std::string_view probabilities_id = "probabilities";
inline constexpr std::string_view annotation_coefficients_id
    = "annotation_coefficients";
inline constexpr std::string_view fitted_values_id = "fitted_values";
inline constexpr std::string_view component_explained_variance_id
    = "component_explained_variance";

template <VarianceLayout Kind>
using marker_variance_dtype_t
    = std::conditional_t<Kind == VarianceLayout::Pooled, double, float>;
template <VarianceLayout Kind>
using marker_variance_writer_t = PayloadWriter<marker_variance_dtype_t<Kind>>;

template <MixtureWeightUpdate Update>
using probability_writer_t = std::conditional_t<
    Update == MixtureWeightUpdate::Enabled,
    PayloadWriter<double>,
    std::monostate>;

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_DRAW_SCHEMA_H_
