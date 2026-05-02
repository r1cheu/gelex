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

#include "gelex/model/bayes/writer/vi_result_writer.h"

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>

#include <fmt/format.h>

#include "gelex/data/reader.h"
#include "gelex/io/text_writer.h"
#include "gelex/model/bayes/writer/covariate_writer.h"
#include "gelex/types/vi_result.h"

namespace gelex
{

vi::ResultWriter::ResultWriter(
    const vi::Result& result,
    const std::filesystem::path& bim_file_path)
    : result_(&result), bim_file_path_(bim_file_path)
{
}

namespace
{

auto write_summary(const vi::Result& result, const std::filesystem::path& path)
    -> void
{
    auto writer = std::make_unique<io::detail::TextWriter>(path);
    writer->write_header({"term", "mean", "stddev"});

    const auto& fixed = result.fixed();
    fixed.for_each_term(
        [&](const std::string& term, Eigen::Index i)
        {
            writer->write(
                fmt::format(
                    "{}\t{}\t{}",
                    term,
                    fixed.coeffs.mean(i),
                    fixed.coeffs.stddev(i)));
        });

    for (const auto& rand : result.random())
    {
        rand.for_each_term(
            [&](const std::string& term, Eigen::Index i)
            {
                writer->write(
                    fmt::format(
                        "{}\t{}\t{}",
                        term,
                        rand.coeffs.mean(i),
                        rand.coeffs.stddev(i)));
            });
    }
}

auto write_snp_effects(
    const vi::Result& result,
    const std::filesystem::path& bim_file_path,
    const std::filesystem::path& path) -> void
{
    const auto* additive = result.genetic(GeneticMode::A);
    if (additive == nullptr)
    {
        return;
    }
    const auto* dominant = result.genetic(GeneticMode::D);

    auto bim = read_bim(bim_file_path);
    auto bim_keys = bim.index().keys();
    auto bim_chrom = bim["chrom"].as<std::string>();
    auto bim_pos = bim["pos"].as<std::int32_t>();
    auto bim_a1 = bim["A1"].as<std::string>();
    auto bim_a2 = bim["A2"].as<std::string>();

    auto writer = std::make_unique<io::detail::TextWriter>(path);

    std::string header = "Index\tID\tChrom\tPosition\tA1\tA2\tA1Freq\tAdd";
    if (dominant != nullptr)
    {
        header += "\tDom";
    }
    writer->write(header);

    for (Eigen::Index i = 0; i < additive->coeffs.size(); ++i)
    {
        auto row = static_cast<std::size_t>(i);
        std::string line = fmt::format("{}\t", i + 1);

        if (i < static_cast<Eigen::Index>(bim.rows()))
        {
            line += fmt::format(
                "{}\t{}\t{}\t{}\t{}",
                bim_keys[row],
                bim_chrom[row],
                bim_pos[row],
                bim_a1[row],
                bim_a2[row]);

            if (result.allele_freq().size() > i)
            {
                line += fmt::format("\t{:.6f}", result.allele_freq()(i));
            }
            else
            {
                line += "\tNA";
            }
        }
        else
        {
            line += "NA\tNA\tNA\tNA\tNA\tNA";
        }

        line += fmt::format("\t{:.6f}", additive->coeffs(i));

        if (dominant != nullptr && i < dominant->coeffs.size())
        {
            line += fmt::format("\t{:.6f}", dominant->coeffs(i));
        }

        writer->write(line);
    }
}

}  // namespace

auto vi::ResultWriter::save(const std::filesystem::path& prefix) const -> void
{
    auto param_path = prefix;
    param_path.replace_extension("param");
    CovariateWriter covar_writer(
        result_->fixed(), result_->random(), param_path);
    covar_writer.write();

    auto summary_path = prefix;
    summary_path.replace_extension("summary");
    write_summary(*result_, summary_path);

    if (!result_->genetics().empty())
    {
        auto snp_path = prefix;
        snp_path.replace_extension(".snp.eff");
        write_snp_effects(*result_, bim_file_path_, snp_path);
    }
}

}  // namespace gelex
