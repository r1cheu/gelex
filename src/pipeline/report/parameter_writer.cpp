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

#include "gelex/pipeline/report/parameter_writer.h"

#include <format>
#include <memory>

#include "gelex/io/text_writer.h"

namespace gelex
{

using Eigen::Index;

ParameterWriter::ParameterWriter(
    const MCMCResult& result,
    const std::filesystem::path& output_path)
    : result_(&result),
      writer_(std::make_unique<detail::TextWriter>(output_path))
{
}

ParameterWriter::~ParameterWriter() = default;

auto ParameterWriter::write() -> void
{
    const auto* additive = result_->additive();
    const auto* dominant = result_->dominant();

    writer_->write_header({"term", "mean", "stddev"});

    write_fixed_effects();
    write_random_effects();
    write_genetic_effect("σ²_add", "h²", additive);
    write_genetic_effect("σ²_dom", "δ²", dominant);
    write_residual_variance();
}

auto ParameterWriter::write_fixed_effects() -> void
{
    const auto* fixed = result_->fixed();
    if (fixed == nullptr)
    {
        return;
    }

    const auto n = fixed->coeffs.size();
    std::vector<std::string> terms(n, "Intercept");
    write_summary_statistics(terms, fixed->coeffs, n);
}

auto ParameterWriter::write_random_effects() -> void
{
    for (const auto& rand : result_->random())
    {
        write_summary_statistics(
            std::vector<std::string>(rand.coeffs.size()),
            rand.coeffs,
            rand.coeffs.size());
        write_summary_statistics(
            std::vector<std::string>(rand.variance.size()),
            rand.variance,
            rand.variance.size());
    }
}

auto ParameterWriter::write_residual_variance() -> void
{
    write_summary_statistics(
        std::vector<std::string>{"σ²_e"},
        result_->residual(),
        result_->residual().size());
}

auto ParameterWriter::write_genetic_effect(
    const std::string& variance_label,
    const std::string& heritability_label,
    const BaseMarkerSummary* effect) -> void
{
    if (effect == nullptr)
    {
        return;
    }

    write_summary_statistics(
        std::vector<std::string>{variance_label},
        effect->variance,
        effect->variance.size());

    std::vector<std::string> proportion_terms;
    proportion_terms.reserve(effect->mixture_proportion.size());
    for (Index i = 0; i < effect->mixture_proportion.size(); ++i)
    {
        proportion_terms.emplace_back(std::format("π[{}]", i));
    }

    write_summary_statistics(
        std::vector<std::string>{heritability_label},
        effect->heritability,
        effect->heritability.size());

    write_summary_statistics(
        proportion_terms,
        effect->mixture_proportion,
        effect->mixture_proportion.size());
}

auto ParameterWriter::write_summary_statistics(
    std::span<const std::string> terms,
    const PosteriorSummary& stats,
    Index n_params) -> void
{
    for (Index i = 0; i < n_params; ++i)
    {
        writer_->write(
            std::format(
                "{}\t{}\t{}", terms[i], stats.mean(i), stats.stddev(i)));
    }
}

}  // namespace gelex
