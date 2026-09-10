// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_SERIALIZATION_IDS_H_
#define GELEX_BAYES_SERIALIZATION_IDS_H_

#include <fmt/compile.h>
#include <fmt/format.h>
#include <string>
#include <string_view>

#include "gelex/genetic_mode.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

inline constexpr std::string_view fixed_coefficients_id = "fixed/coefficients";
inline constexpr std::string_view residual_variance_id = "residual/variance";

inline auto random_coefficients_id(std::string_view name) -> std::string
{
    return fmt::format("random/{}/coefficients", name);
}

inline auto random_variance_id(std::string_view name) -> std::string
{
    return fmt::format("random/{}/variance", name);
}

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

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_SERIALIZATION_IDS_H_
