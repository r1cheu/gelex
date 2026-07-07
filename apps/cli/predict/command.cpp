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

#include <cstddef>
#include <filesystem>
#include <optional>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/data/bed.h"
#include "gelex/data/locus_encoding.h"
#include "gelex/data/reader.h"
#include "gelex/data/snp_alignment.h"
#include "gelex/exception.h"
#include "gelex/io/snpstats.h"

#include "compute.h"
#include "io.h"
#include "reporter.h"
#include "types.h"

auto predict_execute(const cli::PredictConfig& config) -> int
{
    cli::PredictReporter reporter;

    const auto& bfile_prefix = config.bfile;
    const auto& gfile_prefix = config.gfile;

    auto snp_effects = gelex::read_snp_effects(gfile_prefix + ".snpeff");
    auto snpstats = gelex::load_snp_stats(gfile_prefix + ".snpstats");

    if (snpstats.empty())
    {
        throw gelex::GelexException(
            ".snpstats file contains neither additive nor dominance stats.");
    }

    cli::SnpEffects effects;
    for (const auto mode : std::views::keys(snpstats))
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
    auto coefficients = cli::read_coefficients(gfile_prefix + ".param");

    auto fam_df = gelex::read_fam(bfile_prefix + ".fam");
    auto bim_df = gelex::read_bim(bfile_prefix + ".bim");
    auto qcovar_path
        = config.qcovar
              ? std::make_optional<std::filesystem::path>(*config.qcovar)
              : std::nullopt;
    auto dcovar_path
        = config.dcovar
              ? std::make_optional<std::filesystem::path>(*config.dcovar)
              : std::nullopt;
    auto covariates
        = cli::read_covariates(qcovar_path, dcovar_path, coefficients, fam_df);

    auto alignment = gelex::build_snp_alignment(snp_effects, bim_df);
    const auto n_snps = static_cast<std::size_t>(snp_effects.rows());
    reporter.show_snp_selection(
        static_cast<std::size_t>(alignment.num_same + alignment.num_flip),
        static_cast<std::size_t>(alignment.num_absent),
        static_cast<std::size_t>(alignment.num_incompatible),
        n_snps,
        bfile_prefix,
        gfile_prefix + ".snpeff");

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

    auto bed = gelex::open_bed(bfile_prefix, fam_df.index());
    const auto dosage = gelex::load_aligned_genotypes(bed, alignment);

    cli::GenotypeData geno;
    for (const auto& [mode, stats] : snpstats)
    {
        Eigen::MatrixXd encoded = dosage;
        gelex::transform_inplace<double>(
            encoded, gelex::build_loci_encoding(stats));
        geno.emplace(mode, std::move(encoded));
    }

    reporter.show_data_loaded(
        static_cast<std::size_t>(fam_df.rows()),
        static_cast<std::size_t>(snp_effects.rows()),
        coefficients.names.size(),
        snpstats.begin()->second.method);

    auto gebv = cli::compute_gebv(geno, effects);
    auto covar = cli::compute_covariate_effects(covariates, coefficients);

    auto sample_keys = fam_df.index().keys();
    std::vector<std::string> sample_ids(sample_keys.begin(), sample_keys.end());

    cli::PredictResult result{
        .sample_ids = std::move(sample_ids),
        .predictions = gebv.total + covar.total,
        .snp_predictions = std::move(gebv.total),
        .snp_components = std::move(gebv.components),
        .covar_predictions = std::move(covar.per_covariate),
        .covar_names = std::move(covar.covar_names)};

    cli::write_predictions(config.out, result);

    reporter.show_results_written(config.out, result.sample_ids.size());

    return 0;
}
