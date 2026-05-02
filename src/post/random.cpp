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

#include <string_view>

#include <fmt/format.h>

#include "gelex/infra/logging/post_event.h"
#include "gelex/io/detail/binary_reader.h"
#include "gelex/post/detail/utils.h"
#include "gelex/types/genetic_effect_type.h"

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
    auto paths = ref.section_paths();
    auto random_prefix = fmt::format("{}/coeff/", EffectType::random());

    std::vector<std::string_view> effect_names;
    for (auto path : paths)
    {
        if (path.starts_with(random_prefix))
        {
            effect_names.push_back(path.substr(random_prefix.size()));
        }
    }

    std::vector<ParameterDiag> diags;
    for (auto name : effect_names)
    {
        auto levels_path
            = fmt::format("{}/levels/{}", EffectType::random(), name);
        auto levels = ref.to_strings(levels_path);

        std::vector<std::string> coeff_names;
        coeff_names.reserve(levels.size());
        for (auto level : levels)
        {
            coeff_names.emplace_back(fmt::format("{}_{}", name, level));
        }

        auto coeff_path
            = fmt::format("{}/coeff/{}", EffectType::random(), name);
        diags.append_range(
            post::detail::summarize_section(
                readers_,
                coeff_path,
                hdpi_threshold_,
                "Parameter",
                coeff_names));

        auto var_path
            = fmt::format("{}/variance/{}", EffectType::random(), name);
        diags.push_back(
            post::detail::summarize_section(
                readers_,
                var_path,
                hdpi_threshold_,
                "Parameter",
                fmt::format("σ²_{}", name)));
    }
    return diags;
}

}  // namespace gelex
