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
#include <fmt/ranges.h>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cli/formatter.h"
#include "cli/report_printer.h"
#include "gelex/data/dataframe/encode.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/snp_alignment.h"

namespace cli
{

auto PredictReporter::show_snp_selection(
    const gelex::AlignmentPlan& alignment,
    std::string_view bfile_path,
    std::string_view snp_effect_path) const -> void
{
    const auto num_matched
        = static_cast<size_t>(alignment.num_same + alignment.num_flip);
    const auto num_missing = static_cast<size_t>(alignment.num_absent);
    const auto num_mismatched = static_cast<size_t>(alignment.num_incompatible);
    const auto num_total = static_cast<size_t>(alignment.train_count);

    cli::printer().block(gelex::section("SNP Alignment:"));
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

auto PredictReporter::show_covariate_level_mismatches(
    const std::vector<std::pair<std::string, gelex::LevelMismatch>>& mismatches)
    const -> void
{
    for (const auto& [col_name, mismatch] : mismatches)
    {
        if (!mismatch.missing_in_levels.empty())
        {
            cli::printer().warn(
                "Covariate '{}': level(s) [{}] in the data have no fitted "
                "coefficient; affected samples are treated as the reference "
                "level (zero contribution).",
                col_name,
                fmt::join(mismatch.missing_in_levels, ", "));
        }
        if (!mismatch.missing_in_data.empty())
        {
            cli::printer().warn(
                "Covariate '{}': fitted level(s) [{}] are absent from the "
                "data; their columns are zero for all samples.",
                col_name,
                fmt::join(mismatch.missing_in_data, ", "));
        }
    }
}

auto PredictReporter::show_data_loaded(
    size_t num_samples,
    size_t num_snps,
    size_t num_covar_terms,
    gelex::GenotypeMethod geno_method) const -> void
{
    cli::printer().block(gelex::section("Dataset Summary:"));
    cli::printer().line("   {:<13}: {} samples", "Samples", num_samples);
    cli::printer().line("   {:<13}: {} markers", "SNPs", num_snps);
    cli::printer().line("   {:<13}: {}", "Covariates", num_covar_terms);
    cli::printer().line(
        "   {:<13}: {}", "Geno method", fmt::format("{}", geno_method));
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
