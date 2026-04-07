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

#include "gelex/post/genetic_posterior_processor.h"

#include <string>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/infra/logging/post_event.h"
#include "gelex/io/binary_reader.h"
#include "gelex/post/posterior_utils.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

GeneticPosteriorProcessor::GeneticPosteriorProcessor(
    std::span<const detail::BinaryReader> readers,
    double hdpi_threshold)
    : readers_{readers}, hdpi_threshold_{hdpi_threshold}
{
}

auto GeneticPosteriorProcessor::process() -> std::vector<ParameterDiag>
{
    const auto& ref = readers_.front();

    std::vector<ParameterDiag> diags;

    for (auto kind : kAllGeneticKinds)
    {
        auto sect = EffectType::from_genetic(kind);
        if (!ref.contains(fmt::format("{}/coeff", sect)))
        {
            continue;
        }

        auto section = fmt::format("{}", kind);

        // group/proportion
        auto prop_path = fmt::format("{}/group/proportion", sect);
        if (ref.contains(prop_path))
        {
            auto n_pi = ref.to_mat<double>(prop_path).rows();
            std::vector<std::string> pi_names;
            pi_names.reserve(static_cast<size_t>(n_pi));
            for (Eigen::Index i{0}; i < n_pi; ++i)
            {
                pi_names.emplace_back(fmt::format("π[{}]", i));
            }
            diags.append_range(
                detail::summarize_section(
                    readers_, prop_path, hdpi_threshold_, section, pi_names));
        }

        // sign/proportion
        auto sign_path = fmt::format("{}/sign/proportion", sect);
        if (ref.contains(sign_path))
        {
            diags.push_back(
                detail::summarize_section(
                    readers_, sign_path, hdpi_threshold_, section, "p(+)"));
        }
    }

    return diags;
}

}  // namespace gelex
