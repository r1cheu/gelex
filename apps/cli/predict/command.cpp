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
#include "gelex/data/locus_encoding.h"
#include "gelex/data/locus_encoding_types.h"
#include "gelex/data/reader.h"
#include "gelex/data/snp_stats.h"
#include "gelex/exception.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/predict_reader.h"
#include "gelex/io/predict_writer.h"
#include "gelex/io/snpstats.h"
#include "gelex/predict/compute.h"
#include "gelex/predict/snp_alignment.h"
#include "gelex/predict/types.h"
#include "gelex/types/genetic_mode.h"
#include "reporter.h"

namespace
{

auto load_snpstats(const std::filesystem::path& path) -> gelex::SnpStatsData
{
    gelex::BinaryReader reader(path.string());
    gelex::SnpStatsData data;
    if (gelex::has_snp_stats(reader, gelex::GeneticMode::A))
    {
        data.add = gelex::read_snp_stats(reader, gelex::GeneticMode::A);
    }
    if (gelex::has_snp_stats(reader, gelex::GeneticMode::D))
    {
        data.dom = gelex::read_snp_stats(reader, gelex::GeneticMode::D);
    }
    return data;
}

auto build_loci_encoding(const gelex::SnpStats& stats) -> gelex::LociEncoding
{
    gelex::LociEncoding encoding;
    const Eigen::Index n_snps = stats.code.cols();
    encoding.loci.reserve(static_cast<std::size_t>(n_snps));
    for (Eigen::Index j = 0; j < n_snps; ++j)
    {
        gelex::LocusEncoding locus;
        locus.column_index = j;
        locus.code = {stats.code(0, j), stats.code(1, j), stats.code(2, j)};
        locus.missing_encoded_value = 0.0;
        locus.valid = true;
        encoding.loci.push_back(locus);
    }
    return encoding;
}

}  // namespace

auto predict_execute(const cli::PredictConfig& config) -> int
{
    cli::PredictReporter reporter;

    const auto& bfile_prefix = config.bfile;
    const auto& gfile_prefix = config.gfile;

    auto snp_effects = gelex::read_snp_effects(gfile_prefix + ".snpeff");
    auto snpstats = load_snpstats(gfile_prefix + ".snpstats");

    const bool enable_add = snpstats.add.has_value();
    const bool enable_dom = snpstats.dom.has_value();
    if (!enable_add && !enable_dom)
    {
        throw gelex::GelexException(
            ".snpstats file contains neither additive nor dominance stats.");
    }
    if (enable_add && !snp_effects.contains("BETA_A"))
    {
        throw gelex::GelexException(
            ".snpstats file contains additive stats, but SNP effects file "
            "does not have 'BETA_A' column.");
    }
    if (enable_dom && !snp_effects.contains("BETA_D"))
    {
        throw gelex::GelexException(
            ".snpstats file contains dominance stats, but SNP effects file "
            "does not have 'BETA_D' column.");
    }

    auto add_effects = enable_add ? std::make_optional<Eigen::VectorXd>(
                                        snp_effects["BETA_A"].to_map<double>())
                                  : std::nullopt;
    auto dom_effects = enable_dom ? std::make_optional<Eigen::VectorXd>(
                                        snp_effects["BETA_D"].to_map<double>())
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
    auto genotype = gelex::load_aligned_genotypes(bed, alignment);

    gelex::GenotypeData geno;
    if (enable_add && enable_dom)
    {
        geno.dom = genotype;
        geno.add = std::move(genotype);
    }
    else if (enable_add)
    {
        geno.add = std::move(genotype);
    }
    else
    {
        geno.dom = std::move(genotype);
    }

    reporter.show_data_loaded(
        static_cast<std::size_t>(fam_df.rows()),
        static_cast<std::size_t>(snp_effects.rows()),
        coefficients.names.size(),
        (enable_add ? snpstats.add : snpstats.dom)->method);

    if (geno.add)
    {
        const auto encoding = build_loci_encoding(*snpstats.add);
        gelex::transform_inplace<double>(*geno.add, encoding);
    }
    if (geno.dom)
    {
        const auto encoding = build_loci_encoding(*snpstats.dom);
        gelex::transform_inplace<double>(*geno.dom, encoding);
    }

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
