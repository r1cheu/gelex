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

#include "command.h"

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <optional>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/data/snp_alignment.h"
#include "gelex/data/snp_lut_io.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "cli/formatter.h"
#include "cli/report_printer.h"
#include "cli/summary.h"
#include "compute.h"
#include "io.h"

auto predict_execute(const cli::PredictConfig& config) -> int
{
    const auto& bfile_prefix = config.bfile;
    const auto& gfile_prefix = config.gfile;

    // Load and parse the trained model.
    auto snp_effects = gelex::read_snp_effects(gfile_prefix + ".snpeff");
    auto snp_luts = gelex::load_snp_luts(gfile_prefix + ".snplut");
    auto param = gelex::read_param(gfile_prefix + ".param");

    if (snp_luts.empty())
    {
        throw gelex::GelexException(
            ".snplut file contains neither additive nor dominance LUTs.");
    }

    gelex::ModeMap<Eigen::VectorXd> effects;
    for (const auto& [mode, luts] : snp_luts)
    {
        const auto column = fmt::format("BETA_{}", mode);
        if (!snp_effects.contains(column))
        {
            throw gelex::GelexException(
                fmt::format(
                    ".snplut file contains a {} LUT, but SNP effects file "
                    "does not have '{}' column.",
                    mode,
                    column));
        }
        if (luts.cols() != snp_effects.rows())
        {
            throw gelex::GelexException(
                fmt::format(
                    ".snplut mode {} has {} SNPs, but SNP effects file has {}",
                    mode,
                    luts.cols(),
                    snp_effects.rows()));
        }
        effects.emplace(mode, snp_effects[column].to_mat<double>());
    }
    auto term_names = param.index().keys();
    const Eigen::VectorXd coefficients = param["mean"].to_map<double>();

    // Load target genotypes and covariates, then intersect to common samples.
    auto bed = gelex::open_bed(bfile_prefix);

    std::optional<gelex::DataFrame<std::string>> qcovar_df;
    std::optional<gelex::DataFrame<std::string>> dcovar_df;
    std::vector<const gelex::DataFrameIndex<std::string>*> indices{
        &bed.sample_index()};
    if (config.qcovar)
    {
        qcovar_df = gelex::read_qcovar(*config.qcovar);
        indices.push_back(&qcovar_df->index());
    }
    if (config.dcovar)
    {
        dcovar_df = gelex::read_dcovar(*config.dcovar);
        indices.push_back(&dcovar_df->index());
    }

    auto common_index = gelex::intersect<std::string>(indices);
    bed.gather(common_index);
    if (qcovar_df)
    {
        qcovar_df->gather(common_index);
    }
    if (dcovar_df)
    {
        dcovar_df->gather(common_index);
    }

    auto [covariates, level_mismatches] = cli::build_covariate_design(
        term_names,
        qcovar_df,
        dcovar_df,
        static_cast<Eigen::Index>(common_index.size()));
    for (const auto& [column_name, mismatch] : level_mismatches)
    {
        if (!mismatch.missing_in_levels.empty())
        {
            cli::printer().warn(
                "Covariate '{}': level(s) [{}] in the data have no fitted "
                "coefficient; affected samples are treated as the reference "
                "level (zero contribution).",
                column_name,
                fmt::join(mismatch.missing_in_levels, ", "));
        }
        if (!mismatch.missing_in_data.empty())
        {
            cli::printer().warn(
                "Covariate '{}': fitted level(s) [{}] are absent from the "
                "data; their columns are zero for all samples.",
                column_name,
                fmt::join(mismatch.missing_in_data, ", "));
        }
    }

    const auto model_snps = static_cast<std::size_t>(snp_effects.rows());
    cli::Summary{"Dataset Summary"}
        .field("Samples", "{}", common_index.size())
        .field("Variants", "{}", bed.num_snps())
        .show();
    cli::Summary{"Model Summary"}
        .field("Fixed terms", "{}", term_names.size())
        .field("Variants", "{}", model_snps)
        .show();

    // Align SNPs to the model, then expand each mode with its training codes.
    auto alignment = gelex::build_snp_alignment(snp_effects, bed.bim());
    const auto matched_snps
        = static_cast<std::size_t>(alignment.num_same + alignment.num_flip);
    cli::Summary{"SNP Alignment"}
        .field(
            "Matched",
            "{}/{}",
            matched_snps,
            static_cast<std::size_t>(alignment.train_count))
        .field("Missing", "{}", alignment.num_absent)
        .field("Mismatched", "{}", alignment.num_incompatible)
        .show();

    const double missing_ratio
        = static_cast<double>(alignment.missing_pos.size())
          / static_cast<double>(model_snps);
    if (missing_ratio > gelex::max_snp_missing_ratio)
    {
        throw gelex::GelexException(
            fmt::format(
                "{}.snpeff too poorly matched to {}.bim: {} of {} SNPs "
                "missing ({:.1f}%), exceeds {:.0f}% limit",
                gfile_prefix,
                bfile_prefix,
                alignment.missing_pos.size(),
                model_snps,
                missing_ratio * 100.0,
                gelex::max_snp_missing_ratio * 100.0));
    }

    gelex::ModeMap<Eigen::MatrixXd> geno;
    for (const auto& [mode, luts] : snp_luts)
    {
        geno.emplace(
            mode, gelex::expand_aligned_genotypes(bed, alignment, luts));
    }

    // Compute predictions.
    auto gebvs = cli::compute_gebv(geno, effects);
    auto covar
        = cli::compute_covariate_effects(covariates, term_names, coefficients);

    Eigen::VectorXd prediction = covariates * coefficients;
    for (const auto& gebv : std::views::values(gebvs))
    {
        prediction += gebv;
    }

    // Write results.
    auto sample_ids = common_index.keys();
    cli::write_predictions(
        config.out + ".pred.tsv", sample_ids, prediction, covar, gebvs);

    cli::printer().block(cli::results_saved(config.out, ".pred.tsv, .log"));

    return 0;
}
