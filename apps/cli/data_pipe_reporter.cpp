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
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/logging/progress_bar.h"

namespace gelex::cli
{

DataPipeReporter::DataPipeReporter()
    : logger_(gelex::logging::get()), progress_info_(create_progress_info())
{
}

auto DataPipeReporter::on_event(const DataPipeSectionEvent& /*event*/) const
    -> void
{
    logger_->info(gelex::section("[Dataset Summary]"));
}

auto DataPipeReporter::on_event(const PhenotypeLoadedEvent& event) const -> void
{
    logger_->info(
        "   Phenotypes : {} samples ('{}')",
        event.pheno_samples,
        event.trait_name);
}

auto DataPipeReporter::on_event(const CovariatesLoadedEvent& event) const
    -> void
{
    std::string parts;
    if (event.num_quantitative_covariates)
    {
        parts += fmt::format(
            "{} quantitative ({})",
            *event.num_quantitative_covariates,
            gelex::format_names(event.quantitative_names));
    }
    if (event.num_quantitative_covariates && event.num_discrete_covariates)
    {
        parts += ", ";
    }
    if (event.num_discrete_covariates)
    {
        parts += fmt::format(
            "{} discrete ({})",
            *event.num_discrete_covariates,
            gelex::format_names(event.discrete_names));
    }
    logger_->info("   Covariates : {}", parts);
}

auto DataPipeReporter::on_event(const IntersectionEvent& event) const -> void
{
    logger_->info("   Intersection : {} common samples", event.common_samples);
}

auto DataPipeReporter::on_event(const GenotypeLoadedEvent& event) const -> void
{
    const auto effective_snps = event.num_snps - event.monomorphic_snps;
    const std::string label = event.is_dominance ? "Dominance" : "Additive";
    const std::string msg = fmt::format(
        "   {:<13}: {} SNPs ({} monomorphic excluded)",
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

auto DataPipeReporter::on_event(const GrmLoadedEvent& event) const -> void
{
    logger_->info(
        "   GRM        : {} samples ({})", event.num_samples, event.type);
}

auto DataPipeReporter::on_event(const GenotypeProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        progress_info_ = create_progress_info();
        progress_info_.display->show();
    }

    progress_info_.progress_info->message(
        fmt::format(
            "  {}/{} SNPs",
            gelex::AbbrNumber(event.current),
            gelex::AbbrNumber(event.total)));

    if (event.done)
    {
        progress_info_.display->done();
        init_progress_ = false;
    }
}

}  // namespace gelex::cli
