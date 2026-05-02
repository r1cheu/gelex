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

#ifndef GELEX_POST_DETAIL_UTILS_H_
#define GELEX_POST_DETAIL_UTILS_H_

#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/infra/logging/post_event.h"
#include "gelex/infra/stats/diagnostics.h"
#include "gelex/io/detail/binary_reader.h"

namespace gelex::post::detail
{

[[nodiscard]] auto assemble_chains(
    std::span<const ::gelex::io::detail::BinaryReader> readers,
    std::string_view section_path) -> stats::Chains;

[[nodiscard]] auto compute_posterior_summaries(
    const stats::Chains& samples,
    double hdpi_threshold) -> std::vector<ParameterDiag>;

[[nodiscard]] auto summarize_section(
    std::span<const ::gelex::io::detail::BinaryReader> readers,
    std::string_view section_path,
    double hdpi_threshold,
    std::string_view section,
    std::string_view name) -> ParameterDiag;

[[nodiscard]] auto summarize_section(
    std::span<const ::gelex::io::detail::BinaryReader> readers,
    std::string_view section_path,
    double hdpi_threshold,
    std::string_view section,
    std::span<const std::string> names) -> std::vector<ParameterDiag>;

}  // namespace gelex::post::detail

#endif  // GELEX_POST_DETAIL_UTILS_H_
