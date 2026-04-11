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

#include "gelex/pipeline/posterior_analysis_engine.h"

#include <ranges>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/data/genotype/genotype_mmap.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/logging/post_event.h"
#include "gelex/io/binary_reader.h"
#include "gelex/model/bayes/genotype_storage.h"
#include "gelex/post/fixed_posterior_processor.h"
#include "gelex/post/genetic_posterior_processor.h"
#include "gelex/post/genetic_variance_posterior_processor.h"
#include "gelex/post/heritability_posterior_processor.h"
#include "gelex/post/random_posterior_processor.h"
#include "gelex/post/residual_posterior_processor.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

PosteriorAnalysisEngine::PosteriorAnalysisEngine(
    std::span<const std::string_view> samples_prefix,
    double hdpi_threshold,
    std::optional<std::string> gfile)
    : hdpi_threshold_{hdpi_threshold}, gfile_{std::move(gfile)}
{
    for (auto prefix : samples_prefix)
    {
        readers_.emplace_back(fmt::format("{}.samples", prefix));
    }
}

auto PosteriorAnalysisEngine::run(const PostObserver& observer) -> void
{
    if (!check_consistency())
    {
        throw GelexException("Inconsistent section paths across sample files");
    }

    std::span<const detail::BinaryReader> readers{readers_};

    auto fixed_diags
        = FixedPosteriorProcessor{readers, hdpi_threshold_}.process();
    auto random_diags
        = RandomPosteriorProcessor{readers, hdpi_threshold_}.process();
    auto genetic_diags
        = GeneticPosteriorProcessor{readers, hdpi_threshold_}.process();
    auto gebv_diags = process_gebv_variance();
    auto residual_diags
        = ResidualPosteriorProcessor{readers, hdpi_threshold_}.process();

    std::vector<ParameterDiag> diags;

    diags.append_range(std::move(fixed_diags));
    diags.append_range(std::move(random_diags));
    diags.append_range(std::move(gebv_diags));
    diags.append_range(std::move(genetic_diags));
    diags.append_range(std::move(residual_diags));

    notify(
        observer,
        DiagnosticsReadyEvent{
            .diags = std::move(diags),
            .n_chains = static_cast<Eigen::Index>(readers_.size()),
            .n_records = readers_.front()
                             .to_map<double>(fmt::format(
                                 "{}/variance", EffectType::residual()))
                             .cols(),
            .hdpi_prob = hdpi_threshold_});
}

auto PosteriorAnalysisEngine::check_consistency() const -> bool
{
    if (readers_.size() <= 1)
    {
        return !readers_.empty();
    }

    const auto reference = readers_.front().section_paths();
    return std::ranges::all_of(
        readers_ | std::views::drop(1),
        [&](const detail::BinaryReader& reader)
        { return reader.section_paths() == reference; });
}

auto PosteriorAnalysisEngine::process_gebv_variance()
    -> std::vector<ParameterDiag>
{
    if (!gfile_)
    {
        return {};
    }

    std::vector<ParameterDiag> diags;

    std::vector<bayes::GenotypeStorage> genotype_storages;
    std::vector<GeneticInput> genetic_inputs;
    genotype_storages.reserve(kAllGeneticModes.size());

    for (auto kind : kAllGeneticModes)
    {
        auto coeff_path
            = fmt::format("{}/coeff", EffectType::from_genetic(kind));
        if (!readers_.front().contains(coeff_path))
        {
            continue;
        }
        auto gbin_path = fmt::format(
            "{}.{}.gbin", *gfile_, genetic_mode::to_file_suffix(kind));

        genotype_storages.emplace_back(GenotypeMap(gbin_path, kind));
        genetic_inputs.push_back({&genotype_storages.back(), kind});
    }

    if (!genetic_inputs.empty())
    {
        auto [gebv_diags, genetic_variances]
            = GeneticVariancePosteriorProcessor{readers_, genetic_inputs, hdpi_threshold_}
                  .process();
        diags.append_range(std::move(gebv_diags));

        std::vector<GeneticMode> active_kinds;
        active_kinds.reserve(genetic_inputs.size());
        for (const auto& input : genetic_inputs)
        {
            active_kinds.push_back(input.kind);
        }
        diags.append_range(
            HeritabilityPosteriorProcessor{
                readers_, hdpi_threshold_, genetic_variances, active_kinds}
                .process());
    }

    return diags;
}

}  // namespace gelex
