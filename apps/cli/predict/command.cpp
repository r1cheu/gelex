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
#include <CLI/CLI.hpp>
#include <Eigen/Core>

#include "gelex/data/bed.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/logging/predict_event.h"
#include "gelex/io/locistats/reader.h"
#include "gelex/io/predict/input_reader.h"
#include "gelex/predict/snp_alignment.h"
#include "gelex/predict/types.h"
#include "gelex/types/genetic_effect_type.h"
#include "io/predict/writer.h"
#include "predict/compute.h"
#include "predict/standardize.h"
#include "reporter.h"

namespace
{

auto load_sbin(const std::filesystem::path& path) -> gelex::predict::SbinData
{
    gelex::LociStatsReader reader(path.string());
    gelex::predict::SbinData data;
    data.add = reader.read(gelex::EffectType::add());
    if (reader.has(gelex::EffectType::dom()))
    {
        data.dom = reader.read(gelex::EffectType::dom());
        data.has_dom = true;
    }
    return data;
}

}  // namespace

auto predict_execute(CLI::App& predict) -> int
{
    cli::PredictReporter reporter;
    auto observer = reporter.as_observer();
    gelex::notify(observer, gelex::PredictBannerEvent{});

    const auto bfile_prefix = predict.get_option("--bfile")->as<std::string>();
    const auto gfile_prefix = predict.get_option("--gfile")->as<std::string>();

    auto snp_effects
        = gelex::predict::read_snp_effects(gfile_prefix + ".snpeff");
    auto sbin = load_sbin(gfile_prefix + ".sbin");

    bool enable_dom{};
    if (sbin.has_dom)
    {
        if (!snp_effects.contains("BETA_D"))
        {
            throw gelex::GelexException(
                "Sbin file contains dominance effects, but SNP effects file "
                "does not have 'BETA_D' column.");
        }
        enable_dom = true;
    }

    Eigen::VectorXd add_effects = snp_effects["BETA_A"].to_map<double>();
    auto dom_effects = enable_dom ? std::make_optional<Eigen::VectorXd>(
                                        snp_effects["BETA_D"].to_mat<double>())
                                  : std::nullopt;
    auto coefficients
        = gelex::predict::read_coefficients(gfile_prefix + ".param");

    gelex::notify(
        observer,
        gelex::PredictParamsLoadedEvent{
            .bfile_prefix = bfile_prefix,
            .gfile_prefix = gfile_prefix,
            .geno_method = sbin.add.method});

    auto fam_df = gelex::read_fam(bfile_prefix + ".fam");
    auto bim_df = gelex::read_bim(bfile_prefix + ".bim");
    auto qcovar_path
        = predict.get_option("--qcovar")->count() > 0
              ? std::make_optional<std::filesystem::path>(
                    predict.get_option("--qcovar")->as<std::string>())
              : std::nullopt;
    auto dcovar_path
        = predict.get_option("--dcovar")->count() > 0
              ? std::make_optional<std::filesystem::path>(
                    predict.get_option("--dcovar")->as<std::string>())
              : std::nullopt;
    auto covariates = gelex::predict::read_covariates(
        qcovar_path, dcovar_path, coefficients, fam_df);

    auto alignment = gelex::predict::build_snp_alignment(snp_effects, bim_df);
    const auto n_snps = static_cast<std::size_t>(snp_effects.rows());
    if (alignment.num_missing > 0 || alignment.num_mismatched > 0)
    {
        gelex::notify(
            observer,
            gelex::PredictSnpSelectionEvent{
                .num_matched
                = n_snps - static_cast<std::size_t>(alignment.num_missing)
                  - static_cast<std::size_t>(alignment.num_mismatched),
                .num_missing = static_cast<std::size_t>(alignment.num_missing),
                .num_mismatched
                = static_cast<std::size_t>(alignment.num_mismatched),
                .num_total = n_snps,
                .bfile_path = bfile_prefix,
                .snp_effect_path = gfile_prefix + ".snpeff"});

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

    gelex::predict::GenotypeData geno;
    if (sbin.has_dom)
    {
        geno.dom = genotype;
    }
    geno.add = std::move(genotype);

    gelex::notify(
        observer,
        gelex::PredictDataLoadedEvent{
            .num_samples = static_cast<std::size_t>(fam_df.rows()),
            .num_snps = static_cast<std::size_t>(snp_effects.rows()),
            .num_covar_terms = coefficients.names.size()});

    gelex::predict::detail::standardize_genotypes(geno, sbin);

    gelex::predict::SnpEffects effects{
        .add = std::move(add_effects), .dom = std::move(dom_effects)};
    auto gebv = gelex::predict::detail::compute_gebv(geno, effects);
    auto covar = gelex::predict::detail::compute_covariate_effects(
        covariates, coefficients);

    auto sample_keys = fam_df.index().keys();
    std::vector<std::string> sample_ids(sample_keys.begin(), sample_keys.end());

    gelex::predict::PredictResult result{
        .sample_ids = std::move(sample_ids),
        .predictions = gebv.total + covar.total,
        .snp_predictions = std::move(gebv.total),
        .add_predictions = std::move(gebv.add_predictions),
        .dom_predictions = std::move(gebv.dom_predictions),
        .covar_predictions = std::move(covar.per_covariate),
        .covar_names = std::move(covar.covar_names)};

    gelex::predict::detail::PredictWriter writer(
        predict.get_option("--out")->as<std::string>());
    writer.write(result);

    gelex::notify(
        observer,
        gelex::PredictResultsWrittenEvent{
            .output_path = predict.get_option("--out")->as<std::string>(),
            .num_samples = result.sample_ids.size()});

    return 0;
}
