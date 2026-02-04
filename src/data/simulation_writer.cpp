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

#include "gelex/data/simulation_writer.h"

#include <format>
#include <fstream>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>

#include <Eigen/Core>

#include "../src/data/parser.h"
#include "../src/utils/formatter.h"
#include "gelex/logger.h"

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

    auto output = detail::open_file<std::ofstream>(output_path, std::ios::out);

    output << "FID\tIID\tphenotype\n";

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(sample_ids.size());
         ++i)
    {
        const std::string_view full_id(sample_ids[i]);
        std::string_view fid = full_id;
        std::string_view iid = full_id;

        if (const auto pos = full_id.find('_'); pos != std::string_view::npos)
        {
            fid = full_id.substr(0, pos);
            iid = full_id.substr(pos + 1);
        }

        output << std::format("{}\t{}\t{}\n", fid, iid, phenotypes[i]);
    }

    auto logger = logging::get();
    logger->info(
        success(" {:<24}: {}", "Phenotypes saved to", output_path.string()));
}

void SimulationWriter::write_causal_effects(
    std::span<const std::string> snp_ids,
    const std::unordered_map<Eigen::Index, CausalEffect>& effects) const
{
    const auto effects_path = causal_path();

    auto output = detail::open_file<std::ofstream>(effects_path, std::ios::out);

    output << "SNP\tadditive_effect\tdominance_effect\tadd_class\tdom_class\n";

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(snp_ids.size()); ++i)
    {
        auto it = effects.find(i);
        if (it == effects.end())
        {
            continue;
        }
        const auto& effect = it->second;
        output << std::format(
            "{}\t{}\t{}\t{}\t{}\n",
            snp_ids[i],
            effect.additive,
            effect.dominance,
            effect.add_class,
            effect.dom_class);
    }

    auto logger = logging::get();
    logger->info(success(
        " {:<24}: {}", "Causal effects saved to", effects_path.string()));
}

}  // namespace gelex
