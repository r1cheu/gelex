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
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/data/bed.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/predict/input_reader.h"
#include "gelex/io/predict/writer.h"
#include "gelex/io/snpstats/reader.h"
#include "gelex/predict/compute.h"
#include "gelex/predict/snp_alignment.h"
#include "gelex/predict/standardize.h"
#include "gelex/predict/types.h"
#include "gelex/types/genetic_mode.h"
#include "reporter.h"

namespace
{

auto load_snpstats(const std::filesystem::path& path) -> gelex::SnpStatsData
{
    gelex::BinaryReader reader(path.string());
    gelex::SnpStatsData data;
    data.add = gelex::read_snp_stats(reader, gelex::GeneticMode::A);
    if (gelex::has_snp_stats(reader, gelex::GeneticMode::D))
    {
        data.dom = gelex::read_snp_stats(reader, gelex::GeneticMode::D);
        data.has_dom = true;
    }
    return data;
}

}  // namespace

auto predict_execute(const cli::PredictConfig& config) -> int
{
    cli::PredictReporter reporter;

    const auto& bfile_prefix = config.bfile;
    const auto& gfile_prefix = config.gfile;

    auto snp_effects = gelex::read_snp_effects(gfile_prefix + ".snpeff");
    auto snpstats = load_snpstats(gfile_prefix + ".snpstats");

    bool enable_dom{};
    if (snpstats.has_dom)
    {
        if (!snp_effects.contains("BETA_D"))
        {
            throw gelex::GelexException(
                ".snpstats file contains dominance effects, but SNP effects "
                "file "
                "does not have 'BETA_D' column.");
        }
        enable_dom = true;
    }

    Eigen::VectorXd add_effects = snp_effects["BETA_A"].to_map<double>();
    auto dom_effects = enable_dom ? std::make_optional<Eigen::VectorXd>(
                                        snp_effects["BETA_D"].to_mat<double>())
                                  : std::nullopt;
    auto coefficients = gelex::read_coefficients(gfile_prefix + ".param");

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
    auto covariates = gelex::read_covariates(
        qcovar_path, dcovar_path, coefficients, fam_df);

    auto alignment = gelex::build_snp_alignment(snp_effects, bim_df);
    const auto n_snps = static_cast<std::size_t>(snp_effects.rows());
    if (alignment.num_missing > 0 || alignment.num_mismatched > 0)
    {
        reporter.show_snp_selection(
            n_snps - static_cast<std::size_t>(alignment.num_missing)
                - static_cast<std::size_t>(alignment.num_mismatched),
            static_cast<std::size_t>(alignment.num_missing),
            static_cast<std::size_t>(alignment.num_mismatched),
            n_snps,
            bfile_prefix,
            gfile_prefix + ".snpeff");

        throw gelex::GelexException(
            fmt::format(
                "{}.snpeff does not match {}.bim: {} missing SNPs, {} "
                "allele mismatches",
                gfile_prefix,
                bfile_prefix,
                alignment.num_missing,
                alignment.num_mismatched));
    }

    auto bed = gelex::open_bed(bfile_prefix, fam_df.index());
    auto genotype = bed.read_snps<double>(alignment.column_map);

    gelex::GenotypeData geno;
    if (snpstats.has_dom)
    {
        geno.dom = genotype;
    }
    geno.add = std::move(genotype);

    reporter.show_data_loaded(
        static_cast<std::size_t>(fam_df.rows()),
        static_cast<std::size_t>(snp_effects.rows()),
        coefficients.names.size(),
        snpstats.add.method);

    gelex::standardize_genotypes(geno, snpstats);

    gelex::SnpEffects effects{
        .add = std::move(add_effects), .dom = std::move(dom_effects)};
    auto gebv = gelex::compute_gebv(geno, effects);
    auto covar = gelex::compute_covariate_effects(covariates, coefficients);

    auto sample_keys = fam_df.index().keys();
    std::vector<std::string> sample_ids(sample_keys.begin(), sample_keys.end());

    gelex::PredictResult result{
        .sample_ids = std::move(sample_ids),
        .predictions = gebv.total + covar.total,
        .snp_predictions = std::move(gebv.total),
        .add_predictions = std::move(gebv.add_predictions),
        .dom_predictions = std::move(gebv.dom_predictions),
        .covar_predictions = std::move(covar.per_covariate),
        .covar_names = std::move(covar.covar_names)};

    gelex::PredictWriter writer(config.out);
    writer.write(result);

    reporter.show_results_written(config.out, result.sample_ids.size());

    return 0;
}
