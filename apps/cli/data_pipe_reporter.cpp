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

#include "data_pipe_reporter.h"

#include <unistd.h>

#include <fmt/format.h>

#include "gelex/infra/logger.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/utils/formatter.h"

namespace gelex::cli
{

DataPipeReporter::DataPipeReporter() : logger_(gelex::logging::get()) {}

auto DataPipeReporter::on_event(const DataPipeSectionEvent& /*event*/) const
    -> void
{
    logger_->info(gelex::section("[Dataset Summary]"));
    logger_->info("");
}

auto DataPipeReporter::on_event(const DataPipePhenoLoadedEvent& event) const
    -> void
{
    logger_->info(
        gelex::success(
            "Phenotypes : {} samples ('{}')",
            event.pheno_samples,
            event.trait_name));
    logger_->info(
        gelex::success("Genotypes  : {} samples", event.geno_samples));
}

auto DataPipeReporter::on_event(const DataPipeCovarsLoadedEvent& event) const
    -> void
{
    std::string parts;
    if (event.qcovar_loaded > 0)
    {
        parts += fmt::format(
            "{} quantitative ({})",
            event.qcovar_loaded,
            gelex::format_names(event.q_names));
    }
    if (event.qcovar_loaded > 0 && event.dcovar_loaded > 0)
    {
        parts += ", ";
    }
    if (event.dcovar_loaded > 0)
    {
        parts += fmt::format(
            "{} discrete ({})",
            event.dcovar_loaded,
            gelex::format_names(event.d_names));
    }
    logger_->info(gelex::success("Covariates : {}", parts));
}

auto DataPipeReporter::on_event(const DataPipeIntersectedEvent& event) const
    -> void
{
    logger_->info(
        gelex::success(
            "Intersection : {} common, {} excluded",
            event.common_samples,
            event.excluded_samples));
}

auto DataPipeReporter::on_event(const DataPipeGenotypeLoadedEvent& event) const
    -> void
{
    const auto effective_snps = event.num_snps - event.monomorphic_snps;
    const std::string label = event.is_dominance ? "Dominance" : "Additive";
    const std::string msg = gelex::success(
        "{:<13}: {} SNPs ({} monomorphic excluded)",
        label,
        gelex::AbbrNumber(effective_snps),
        gelex::AbbrNumber(event.monomorphic_snps));

    if (isatty(fileno(stdout)) != 0)
    {
        logger_->info("{}", "\033[A\r" + msg + "\033[K");
    }
    else
    {
        logger_->info("{}", msg);
    }
}

}  // namespace gelex::cli
