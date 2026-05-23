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

#include "gelex/post/random.h"

#include <string>
#include <string_view>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/infra/logging/post_event.h"
#include "gelex/io/detail/binary_reader.h"
#include "gelex/post/detail/utils.h"

namespace gelex
{

RandomPosteriorProcessor::RandomPosteriorProcessor(
    std::span<const io::detail::BinaryReader> readers,
    double hdpi_threshold)
    : readers_{readers}, hdpi_threshold_{hdpi_threshold}
{
}

auto RandomPosteriorProcessor::process() -> std::vector<ParameterDiag>
{
    const auto& ref = readers_.front();

    std::vector<ParameterDiag> diags;
    for (size_t index = 0;; ++index)
    {
        const auto coeff_path = fmt::format("random/{}/coeffs", index);
        if (!ref.contains(coeff_path))
        {
            break;
        }

        std::vector<std::string> coeff_names;
        const auto n_params = ref.to_map<double>(coeff_path).rows();
        coeff_names.reserve(static_cast<size_t>(n_params));
        for (Eigen::Index i = 0; i < n_params; ++i)
        {
            coeff_names.push_back(fmt::format("random[{}][{}]", index, i));
        }

        diags.append_range(
            post::detail::summarize_section(
                readers_,
                coeff_path,
                hdpi_threshold_,
                "Parameter",
                coeff_names));

        const auto var_path = fmt::format("random/{}/variance", index);
        diags.push_back(
            post::detail::summarize_section(
                readers_,
                var_path,
                hdpi_threshold_,
                "Parameter",
                fmt::format("σ²_random[{}]", index)));
    }
    return diags;
}

}  // namespace gelex
