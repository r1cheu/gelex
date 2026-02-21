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

#include "gelex/io/sim/simulation_writer.h"

#include <algorithm>
#include <format>
#include <span>
#include <string>
#include <string_view>

#include <Eigen/Core>

#include "gelex/io/text_writer.h"
#include "gelex/types/sample_id.h"

namespace gelex
{

SimulationWriter::SimulationWriter(std::filesystem::path output_prefix)
    : output_prefix_(std::move(output_prefix))
{
}

auto SimulationWriter::phenotype_path() const -> std::filesystem::path
{
    auto path = output_prefix_;
    path.replace_extension(".phen");
    return path;
}

auto SimulationWriter::causal_path() const -> std::filesystem::path
{
    auto path = output_prefix_;
    path.replace_extension(".causal");
    return path;
}

void SimulationWriter::write_phenotypes(
    const Eigen::Ref<const Eigen::VectorXd>& phenotypes,
    std::span<const std::string> sample_ids) const
{
    const auto output_path = phenotype_path();

    detail::TextWriter writer(output_path);

    writer.write_header({"FID", "IID", "phenotype"});

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(sample_ids.size());
         ++i)
    {
        const std::string_view full_id(sample_ids[i]);
        auto [fid, iid] = split_sample_id(full_id);
        writer.write(std::format("{}\t{}\t{}", fid, iid, phenotypes[i]));
    }
}

void SimulationWriter::write_causal_effects(
    std::span<const std::string> snp_ids,
    const CausalEffects& effects) const
{
    const auto effects_path = causal_path();

    detail::TextWriter writer(effects_path);

    writer.write_header(
        {"SNP",
         "additive_effect",
         "dominance_effect",
         "add_class",
         "dom_class"});

    const Eigen::Index effect_count = effects.size();
    const Eigen::Index row_count
        = std::min(static_cast<Eigen::Index>(snp_ids.size()), effect_count);

    for (Eigen::Index i = 0; i < row_count; ++i)
    {
        writer.write(
            std::format(
                "{}\t{}\t{}\t{}\t{}",
                snp_ids[i],
                effects.additive(i),
                effects.dominance(i),
                effects.add_class(i),
                effects.dom_class(i)));
    }
}

}  // namespace gelex
