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
#include <cstddef>
#include <string>
#include <string_view>

#include "cli/report_printer.h"
#include "gelex/data/genotype_method.h"
#include "gelex/infra/logging/formatter.h"
#include "version.h"

namespace cli
{

auto PredictReporter::show_banner() const -> void
{
    cli::printer().block(
        gelex::command_banner(PROJECT_VERSION, "Genomic Prediction"));
}

auto PredictReporter::show_params_loaded(
    std::string_view bfile_prefix,
    std::string_view gfile_prefix,
    gelex::GenotypeMethod geno_method) const -> void
{
    cli::printer().block(gelex::section("[Config]"));
    cli::printer().line("  {:<12}: {}", "bfile", bfile_prefix);
    cli::printer().line("  {:<12}: {}", "gfile", gfile_prefix);
    cli::printer().line(
        "  {:<12}: {}", "Geno method", fmt::format("{}", geno_method));
}

auto PredictReporter::show_snp_selection(
    size_t num_matched,
    size_t num_missing,
    size_t num_mismatched,
    size_t num_total,
    std::string_view bfile_path,
    std::string_view snp_effect_path) const -> void
{
    cli::printer().block(gelex::section("[SNP Alignment]"));
    cli::printer().line("   {:<13}: {}/{}", "Matched", num_matched, num_total);
    cli::printer().line("   {:<13}: {}", "Missing", num_missing);
    cli::printer().line("   {:<13}: {}", "Mismatched", num_mismatched);

    if (num_mismatched > 0)
    {
        std::string plink_hint = fmt::format(
            "plink2 --bfile {} --alt1-allele {} 4 1 1 --make-bed --out "
            "<output>",
            bfile_path,
            snp_effect_path);
        cli::printer().warn(
            "Allele mismatch detected for {} SNPs. "
            "To fix, run:\n  {}",
            num_mismatched,
            plink_hint);
    }
}

auto PredictReporter::show_data_loaded(
    size_t num_samples,
    size_t num_snps,
    size_t num_covar_terms) const -> void
{
    cli::printer().block(gelex::section("[Dataset Summary]"));
    cli::printer().line("   {:<13}: {} samples", "Samples", num_samples);
    cli::printer().line("   {:<13}: {} markers", "SNPs", num_snps);
    cli::printer().line("   {:<13}: {}", "Covariates", num_covar_terms);
}

auto PredictReporter::show_results_written(
    std::string_view output_path,
    size_t num_samples) const -> void
{
    cli::printer().block(
        gelex::success(
            "Results saved to '{}' ({} samples)", output_path, num_samples));
}

}  // namespace cli
