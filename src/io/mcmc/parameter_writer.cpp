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

#include "gelex/io/mcmc/parameter_writer.h"

#include <fmt/format.h>
#include <memory>
#include <string>

#include "gelex/io/detail/text_writer.h"
#include "gelex/model/bayes/labels.h"

namespace gelex
{

using Eigen::Index;

ParameterWriter::ParameterWriter(
    const mcmc::Result& result,
    const std::filesystem::path& output_path)
    : result_(&result),
      writer_(std::make_unique<io::detail::TextWriter>(output_path))
{
}

ParameterWriter::~ParameterWriter() = default;

auto ParameterWriter::write() -> void
{
    writer_->write_header({"term", "mean", "stddev"});

    write_fixed_effects();
    write_random_effects();
    for (const auto& summary : result_->genetics())
    {
        write_genetic_effect(summary);
    }
    write_residual_variance();
}

auto ParameterWriter::write_fixed_effects() -> void
{
    const auto& fixed = result_->fixed();
    fixed.for_each_term(
        [&](const std::string& term, Eigen::Index i)
        {
            writer_->write(
                fmt::format(
                    "{}\t{}\t{}",
                    term,
                    fixed.coeffs.mean(i),
                    fixed.coeffs.stddev(i)));
        });
}

auto ParameterWriter::write_random_effects() -> void
{
    for (const auto& rand : result_->random())
    {
        rand.for_each_term(
            [&](const std::string& term, Eigen::Index i)
            {
                writer_->write(
                    fmt::format(
                        "{}\t{}\t{}",
                        term,
                        rand.coeffs.mean(i),
                        rand.coeffs.stddev(i)));
            });

        writer_->write(
            fmt::format(
                "σ²_{}\t{}\t{}",
                rand.name,
                rand.variance.mean(0),
                rand.variance.stddev(0)));
    }
}

auto ParameterWriter::write_residual_variance() -> void
{
    write_summary_statistics(
        std::vector<std::string>{"σ²_e"}, result_->residual());
}

auto ParameterWriter::write_genetic_effect(const GeneticSummary& effect) -> void
{
    write_summary_statistics(
        std::vector<std::string>{
            std::string{bayes::to_variance_label(effect.type)}},
        effect.variance);

    write_summary_statistics(
        std::vector<std::string>{
            std::string{bayes::to_heritability_label(effect.type)}},
        effect.heritability);

    if (effect.group)
    {
        const auto& base = assignment(*effect.group);
        std::vector<std::string> proportion_terms;
        proportion_terms.reserve(base.mixture_proportion.size());
        for (Index i = 0; i < base.mixture_proportion.size(); ++i)
        {
            proportion_terms.emplace_back(fmt::format("π[{}]", i));
        }

        write_summary_statistics(proportion_terms, base.mixture_proportion);
    }
}

auto ParameterWriter::write_summary_statistics(
    std::span<const std::string> terms,
    const PosteriorSummary& stats) -> void
{
    for (Index i = 0; i < stats.size(); ++i)
    {
        writer_->write(
            fmt::format(
                "{}\t{}\t{}", terms[i], stats.mean(i), stats.stddev(i)));
    }
}

}  // namespace gelex
