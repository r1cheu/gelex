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

#include <fmt/format.h>

#include "gelex/infra/logging/post_event.h"
#include "gelex/io/detail/binary_reader.h"
#include "gelex/post/detail/utils.h"
#include "gelex/types/genetic_effect_type.h"

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
    auto names_path = fmt::format("{}/names", EffectType::fixed());
    if (!ref.contains(names_path))
    {
        return {};
    }

    auto names = ref.to_strings(names_path);

    std::vector<std::string> row_names;
    for (auto name : names)
    {
        auto levels_path
            = fmt::format("{}/levels/{}", EffectType::fixed(), name);
        auto levels = ref.to_strings(levels_path);
        if (levels.size() == 1 && levels[0] == name)
        {
            row_names.emplace_back(name);
        }
        else
        {
            for (auto level : levels)
            {
                row_names.emplace_back(fmt::format("{}_{}", name, level));
            }
        }
    }

    auto coeff_path = fmt::format("{}/coeff", EffectType::fixed());
    return post::detail::summarize_section(
        readers_, coeff_path, hdpi_threshold_, "Parameter", row_names);
}

}  // namespace gelex
