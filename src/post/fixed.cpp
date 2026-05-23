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

#include "gelex/post/fixed.h"

#include <string>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/infra/logging/post_event.h"
#include "gelex/io/detail/binary_reader.h"
#include "gelex/post/detail/utils.h"

namespace gelex
{

FixedPosteriorProcessor::FixedPosteriorProcessor(
    std::span<const io::detail::BinaryReader> readers,
    double hdpi_threshold)
    : readers_{readers}, hdpi_threshold_{hdpi_threshold}
{
}

auto FixedPosteriorProcessor::process() -> std::vector<ParameterDiag>
{
    const auto& ref = readers_.front();
    constexpr std::string_view coeff_path = "fixed/0/coeffs";
    if (!ref.contains(coeff_path))
    {
        return {};
    }

    std::vector<std::string> row_names;
    const auto n_params = ref.to_map<double>(coeff_path).rows();
    row_names.reserve(static_cast<size_t>(n_params));
    for (Eigen::Index i = 0; i < n_params; ++i)
    {
        row_names.push_back(fmt::format("fixed[{}]", i));
    }

    return post::detail::summarize_section(
        readers_, coeff_path, hdpi_threshold_, "Parameter", row_names);
}

}  // namespace gelex
