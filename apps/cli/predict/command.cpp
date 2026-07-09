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
#include <optional>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/locus_encoding.h"
#include "gelex/data/reader.h"
#include "gelex/data/snp_alignment.h"
#include "gelex/exception.h"
#include "gelex/io/snpstats.h"

#include "compute.h"
#include "io.h"
#include "reporter.h"

auto predict_execute(const cli::PredictConfig& config) -> int
{
    cli::PredictReporter reporter;

    const auto& bfile_prefix = config.bfile;
    const auto& gfile_prefix = config.gfile;

    // Load and parse the trained model.
    auto snp_effects = gelex::read_snp_effects(gfile_prefix + ".snpeff");
    auto snp_stats = gelex::load_snp_stats(gfile_prefix + ".snpstats");
    auto param = gelex::read_param(gfile_prefix + ".param");

    if (snp_stats.empty())
    {
        throw gelex::GelexException(
            ".snpstats file contains neither additive nor dominance stats.");
    }

    gelex::ModeMap<Eigen::VectorXd> effects;
    for (const auto mode : std::views::keys(snp_stats))
    {
        const auto column = fmt::format("BETA_{}", mode);
        if (!snp_effects.contains(column))
        {
            throw gelex::GelexException(
                fmt::format(
                    ".snpstats file contains {} stats, but SNP effects file "
                    "does not have '{}' column.",
                    mode,
                    column));
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
    reporter.show_covariate_level_mismatches(level_mismatches);

    // Align SNPs to the model, then load the aligned dosage.
    auto alignment = gelex::build_snp_alignment(snp_effects, bed.bim());
    reporter.show_snp_selection(alignment);

    const auto n_snps = static_cast<std::size_t>(snp_effects.rows());
    const double missing_ratio
        = static_cast<double>(alignment.missing_pos.size())
          / static_cast<double>(n_snps);
    if (missing_ratio > gelex::MAX_SNP_MISSING_RATIO)
    {
        throw gelex::GelexException(
            fmt::format(
                "{}.snpeff too poorly matched to {}.bim: {} of {} SNPs "
                "missing ({:.1f}%), exceeds {:.0f}% limit",
                gfile_prefix,
                bfile_prefix,
                alignment.missing_pos.size(),
                n_snps,
                missing_ratio * 100.0,
                gelex::MAX_SNP_MISSING_RATIO * 100.0));
    }

    const auto dosage = gelex::load_aligned_genotypes(bed, alignment);

    gelex::ModeMap<Eigen::MatrixXd> geno;
    for (const auto& [mode, stats] : snp_stats)
    {
        Eigen::MatrixXd encoded = dosage;
        gelex::transform_inplace<double>(
            encoded, gelex::build_loci_encoding(stats));
        geno.emplace(mode, std::move(encoded));
    }

    reporter.show_data_loaded(
        static_cast<std::size_t>(common_index.size()),
        n_snps,
        term_names.size());

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
    cli::write_predictions(config.out, sample_ids, prediction, covar, gebvs);

    reporter.show_results_written(config.out, sample_ids.size());

    return 0;
}
