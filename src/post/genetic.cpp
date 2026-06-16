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

#include "gelex/post/genetic.h"

#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/io/binary_reader.h"
#include "gelex/post/detail/utils.h"
#include "gelex/post/diagnostic.h"

namespace gelex
{

namespace
{

auto block_section_name(const io::BinaryReader& reader, std::size_t block_index)
    -> std::string
{
    const auto modes_path = fmt::format("genetic_block/{}/modes", block_index);
    const auto modes = reader.to_strings(modes_path);
    if (modes.empty())
    {
        return fmt::format("genetic_block[{}]", block_index);
    }
    if (modes.size() == 1)
    {
        return std::string(modes[0]);
    }

    std::string name;
    for (std::size_t i = 0; i < modes.size(); ++i)
    {
        if (i > 0)
        {
            name += ",";
        }
        name += modes[i];
    }
    return name;
}

}  // namespace

GeneticPosteriorProcessor::GeneticPosteriorProcessor(
    std::span<const io::BinaryReader> readers,
    double hdpi_threshold)
    : readers_{readers}, hdpi_threshold_{hdpi_threshold}
{
}

auto GeneticPosteriorProcessor::process() -> std::vector<ParameterDiag>
{
    const auto& ref = readers_.front();

    std::vector<ParameterDiag> diags;
    for (std::size_t block_index = 0;; ++block_index)
    {
        const auto modes_path
            = fmt::format("genetic_block/{}/modes", block_index);
        if (!ref.contains(modes_path))
        {
            break;
        }

        const auto section = block_section_name(ref, block_index);
        for (std::size_t slot = 0;; ++slot)
        {
            const auto assignment_path = fmt::format(
                "genetic_block/{}/prior_state/proportion/{}/assignment",
                block_index,
                slot);
            const auto prop_path = fmt::format(
                "genetic_block/{}/prior_state/proportion/{}/value",
                block_index,
                slot);
            if (!ref.contains(assignment_path) && !ref.contains(prop_path))
            {
                break;
            }
            if (ref.contains(prop_path))
            {
                const auto n_pi = ref.to_mat<double>(prop_path).rows();
                std::vector<std::string> pi_names;
                pi_names.reserve(static_cast<size_t>(n_pi));
                for (Eigen::Index i{0}; i < n_pi; ++i)
                {
                    pi_names.push_back(fmt::format("π[{}]", i));
                }
                diags.append_range(
                    post::detail::summarize_section(
                        readers_,
                        prop_path,
                        hdpi_threshold_,
                        section,
                        pi_names));
            }
        }
    }

    return diags;
}

}  // namespace gelex
