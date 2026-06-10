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

#include "reporter.h"

#include <fmt/format.h>
#include <string>

#include "cli/report_printer.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/logging/predict_event.h"
#include "version.h"

namespace cli
{

auto PredictReporter::on_event(const gelex::PredictBannerEvent& /*event*/) const
    -> void
{
    cli::printer().block(
        gelex::command_banner(PROJECT_VERSION, "Genomic Prediction"));
}

auto PredictReporter::on_event(
    const gelex::PredictParamsLoadedEvent& event) const -> void
{
    cli::printer().block(gelex::section("[Config]"));
    cli::printer().line("  {:<12}: {}", "bfile", event.bfile_prefix);
    cli::printer().line("  {:<12}: {}", "gfile", event.gfile_prefix);
    cli::printer().line(
        "  {:<12}: {}", "Geno method", fmt::format("{}", event.geno_method));
}

auto PredictReporter::on_event(
    const gelex::PredictSnpSelectionEvent& event) const -> void
{
    cli::printer().block(gelex::section("[SNP Alignment]"));
    cli::printer().line(
        "   {:<13}: {}/{}", "Matched", event.num_matched, event.num_total);
    cli::printer().line("   {:<13}: {}", "Missing", event.num_missing);
    cli::printer().line("   {:<13}: {}", "Mismatched", event.num_mismatched);

    if (event.num_mismatched > 0)
    {
        std::string plink_hint = fmt::format(
            "plink2 --bfile {} --alt1-allele {} 4 1 1 --make-bed --out "
            "<output>",
            event.bfile_path,
            event.snp_effect_path);
        cli::printer().warn(
            "Allele mismatch detected for {} SNPs. "
            "To fix, run:\n  {}",
            event.num_mismatched,
            plink_hint);
    }
}

auto PredictReporter::on_event(const gelex::PredictDataLoadedEvent& event) const
    -> void
{
    cli::printer().block(gelex::section("[Dataset Summary]"));
    cli::printer().line("   {:<13}: {} samples", "Samples", event.num_samples);
    cli::printer().line("   {:<13}: {} markers", "SNPs", event.num_snps);
    cli::printer().line("   {:<13}: {}", "Covariates", event.num_covar_terms);
}

auto PredictReporter::on_event(
    const gelex::PredictResultsWrittenEvent& event) const -> void
{
    cli::printer().block(
        gelex::success(
            "Results saved to '{}' ({} samples)",
            event.output_path,
            event.num_samples));
}

}  // namespace cli
