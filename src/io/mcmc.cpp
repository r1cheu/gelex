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

#include "gelex/io/mcmc.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <ranges>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/algo/mcmc/result.h"
#include "gelex/bayes/model.h"
#include "gelex/data/reader.h"
#include "gelex/infra/stats/result.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex::mcmc
{

auto write_params(const Result& result, std::string_view prefix) -> void
{
    io::detail::TextWriter writer(fmt::format("{}.param", prefix));
    writer.write_header({"term", "mean", "stddev"});

    for (const auto& record : result.records())
    {
        const std::string_view path{record.path};
        if (path != "state/fixed/coeffs"
            && !(
                path.starts_with("state/random_") && path.ends_with("/coeffs")))
        {
            continue;
        }
        if (!record.names)
        {
            continue;
        }
        if (!std::holds_alternative<stats::RunningStatsResult>(record.value))
        {
            continue;
        }

        const auto& stats = std::get<stats::RunningStatsResult>(record.value);
        const auto& names = *record.names;
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(names.size());
             ++i)
        {
            writer.write(
                fmt::format(
                    "{}\t{:.8e}\t{:.8e}",
                    names[static_cast<std::size_t>(i)],
                    stats.mean(i),
                    stats.stddev(i)));
        }
    }
}

auto write_summary(const Result& result, std::string_view prefix) -> void
{
    io::detail::TextWriter writer(fmt::format("{}.summary", prefix));
    writer.write_header({"term", "effect", "mean", "stddev"});

    for (const auto& record : result.records())
    {
        const std::string_view path{record.path};
        if (path.ends_with("/pve"))
        {
            continue;
        }
        if (!record.names)
        {
            continue;
        }

        if (!std::holds_alternative<stats::RunningStatsResult>(record.value))
        {
            continue;
        }

        const auto& stats = std::get<stats::RunningStatsResult>(record.value);
        const auto& names = *record.names;
        std::string_view effect{"-"};
        std::size_t start{};
        while (start < path.size())
        {
            const auto end = path.find('/', start);
            const auto segment = end == std::string_view::npos
                                     ? path.substr(start)
                                     : path.substr(start, end - start);
            if (segment == "A" || segment == "D")
            {
                effect = segment;
                break;
            }
            if (end == std::string_view::npos)
            {
                break;
            }
            start = end + 1;
        }
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(names.size());
             ++i)
        {
            writer.write(
                fmt::format(
                    "{}\t{}\t{:.8e}\t{:.8e}",
                    names[static_cast<std::size_t>(i)],
                    effect,
                    stats.mean(i),
                    stats.stddev(i)));
        }
    }
}

auto write_snp_eff(
    const Result& result,
    const BayesModel& model,
    const std::filesystem::path& bim_path,
    std::string_view prefix) -> void
{
    std::unordered_map<std::string, std::size_t> column_indices;
    std::vector<const Eigen::VectorXd*> column_values;
    auto register_column = [&column_indices, &column_values](
                               std::string column, const Eigen::VectorXd& value)
    {
        const auto it = column_indices.find(column);
        if (it != column_indices.end())
        {
            column_values[it->second] = &value;
            return;
        }

        column_indices.emplace(std::move(column), column_values.size());
        column_values.push_back(&value);
    };

    for (const auto& record : result.records())
    {
        if (!std::holds_alternative<stats::RunningStatsResult>(record.value))
        {
            continue;
        }

        const std::string_view path{record.path};
        const auto& stats = std::get<stats::RunningStatsResult>(record.value);
        if (path == "state/genetic/pve")
        {
            register_column("PVE", stats.mean);
        }
        else if (path.ends_with("/genetic/coeffs") && path.contains("/A/"))
        {
            register_column("BETA_A", stats.mean);
            register_column("SE_A", stats.stddev);
        }
        else if (path.ends_with("/genetic/coeffs") && path.contains("/D/"))
        {
            register_column("BETA_D", stats.mean);
            register_column("SE_D", stats.stddev);
        }
        else if (path.ends_with("/genetic/pve") && path.contains("/A/"))
        {
            register_column("PVE_A", stats.mean);
        }
        else if (path.ends_with("/genetic/pve") && path.contains("/D/"))
        {
            register_column("PVE_D", stats.mean);
        }
        else if (path.ends_with("/pip") && path.contains("/A/"))
        {
            register_column("PIP_A", stats.mean);
        }
        else if (path.ends_with("/pip") && path.contains("/D/"))
        {
            register_column("PIP_D", stats.mean);
        }
    }

    const bool has_additive = column_indices.contains("BETA_A");
    const bool has_dominance = column_indices.contains("BETA_D");
    const bool needs_total_pve = has_additive && has_dominance;

    const auto* allele_design = has_additive ? model.genetic(GeneticMode::A)
                                             : model.genetic(GeneticMode::D);
    const auto& A1freq = allele_design->X.A1freq();
    const auto n_snps = A1freq.size();

    auto bim = read_bim(bim_path);

    std::vector<std::string> columns{"CHR", "SNP", "BP", "A1", "A2", "A1FREQ"};
    if (has_additive)
    {
        columns.insert(columns.end(), {"BETA_A", "SE_A", "PVE_A"});
        if (column_indices.contains("PIP_A"))
        {
            columns.emplace_back("PIP_A");
        }
    }
    if (has_dominance)
    {
        columns.insert(columns.end(), {"BETA_D", "SE_D", "PVE_D"});
        if (column_indices.contains("PIP_D"))
        {
            columns.emplace_back("PIP_D");
        }
    }
    if (needs_total_pve)
    {
        columns.emplace_back("PVE");
    }

    std::string header;
    for (const auto [i, column] : std::views::enumerate(columns))
    {
        if (i != 0)
        {
            header.push_back('\t');
        }
        header += column;
    }

    io::detail::TextWriter writer(fmt::format("{}.snpeff", prefix));
    writer.write(header);

    std::vector<const Eigen::VectorXd*> ordered_values;
    ordered_values.reserve(columns.size() - 6);
    for (const auto& column : columns | std::views::drop(6))
    {
        ordered_values.push_back(column_values[column_indices.at(column)]);
    }

    const auto keys = bim.index().keys();
    const auto chrom = bim["chrom"].as<std::string>();
    const auto pos = bim["pos"].as<std::int32_t>();
    const auto a1 = bim["A1"].as<std::string>();
    const auto a2 = bim["A2"].as<std::string>();
    std::string line;
    line.reserve(128 + (ordered_values.size() * 16));
    for (Eigen::Index i = 0; i < n_snps; ++i)
    {
        const auto row = static_cast<std::size_t>(i);
        line.clear();
        fmt::format_to(
            std::back_inserter(line),
            "{}\t{}\t{}\t{}\t{}\t{:.8e}",
            chrom[row],
            keys[row],
            pos[row],
            a1[row],
            a2[row],
            A1freq(i));

        for (const auto* value : ordered_values)
        {
            fmt::format_to(std::back_inserter(line), "\t{:.8e}", (*value)(i));
        }
        writer.write(line);
    }
}

}  // namespace gelex::mcmc
