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

#include "gelex/pipeline/predict_engine.h"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/genotype/bed_pipe.h"
#include "gelex/data/reader.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/locistats_reader.h"
#include "gelex/predict/input_reader.h"
#include "gelex/predict/snp_alignment.h"
#include "predict/compute.h"
#include "predict/standardize.h"
#include "predict/writer.h"

namespace gelex
{

PredictEngine::PredictEngine(Config config) : config_(std::move(config)) {}

auto PredictEngine::load_sbin(const std::filesystem::path& path)
    -> predict::SbinData
{
    LociStatsReader reader(path.string());
    predict::SbinData data;
    data.add = reader.read(EffectType::add());
    if (reader.has(EffectType::dom()))
    {
        data.dom = reader.read(EffectType::dom());
        data.has_dom = true;
    }
    return data;
}

auto PredictEngine::load_params() const -> PredictParams
{
    auto snp_effects
        = predict::read_snp_effects(config_.gfile_prefix + ".snp.eff");
    auto sbin = load_sbin(config_.gfile_prefix + ".sbin");

    bool enable_dom{};
    if (sbin.has_dom)
    {
        if (!snp_effects.contains("Dom"))
        {
            throw GelexException(
                "Sbin file contains dominance effects, but SNP effects file "
                "does not have 'Dom' column.");
        }
        enable_dom = true;
    }

    Eigen::VectorXd add_effects = snp_effects["Add"].to_map<double>();
    auto dom_effects = enable_dom ? std::make_optional<Eigen::VectorXd>(
                                        snp_effects["Dom"].to_mat<double>())
                                  : std::nullopt;

    return PredictParams{
        .snp_effects = std::move(snp_effects),
        .add_effects = std::move(add_effects),
        .dom_effects = std::move(dom_effects),
        .coefficients
        = predict::read_coefficients(config_.gfile_prefix + ".param"),
        .sbin = std::move(sbin)};
}

auto PredictEngine::load_data(const PredictParams& params) const -> PredictData
{
    auto fam_path = config_.bed_path;
    fam_path.replace_extension(".fam");
    auto bim_path = config_.bed_path;
    bim_path.replace_extension(".bim");

    auto fam_df = read_fam(fam_path);
    auto bim_df = read_bim(bim_path);
    auto covariates = predict::read_covariates(
        config_.qcovar_path, config_.dcovar_path, params.coefficients, fam_df);

    return PredictData{
        .fam_df = std::move(fam_df),
        .bim_df = std::move(bim_df),
        .covariates = std::move(covariates)};
}

auto PredictEngine::align(
    const PredictParams& params,
    const PredictData& data,
    const PredictObserver& observer) const -> predict::SnpAlignment
{
    auto alignment
        = predict::build_snp_alignment(params.snp_effects, data.bim_df);
    const auto n_snps = static_cast<size_t>(params.snp_effects.rows());

    if (alignment.num_missing > 0 || alignment.num_mismatched > 0)
    {
        notify(
            observer,
            PredictSnpSelectionEvent{
                .num_matched = n_snps
                               - static_cast<size_t>(alignment.num_missing)
                               - static_cast<size_t>(alignment.num_mismatched),
                .num_missing = static_cast<size_t>(alignment.num_missing),
                .num_mismatched = static_cast<size_t>(alignment.num_mismatched),
                .num_total = n_snps,
                .bfile_path = config_.bfile_prefix,
                .snp_effect_path = config_.gfile_prefix + ".snp.eff"});
    }

    return alignment;
}

auto PredictEngine::select(
    const PredictData& data,
    const predict::SnpAlignment& alignment,
    bool has_dom) const -> predict::GenotypeData
{
    genotype::BedPipe bed_pipe(config_.bed_path, data.fam_df.index());
    auto genotype = bed_pipe.select(alignment.column_map);

    predict::GenotypeData geno;
    if (has_dom)
    {
        geno.dom = genotype;
    }
    geno.add = std::move(genotype);
    return geno;
}

auto PredictEngine::run(const PredictObserver& observer) -> void
{
    auto params = load_params();

    notify(
        observer,
        PredictParamsLoadedEvent{
            .bfile_prefix = config_.bfile_prefix,
            .gfile_prefix = config_.gfile_prefix,
            .geno_method = params.sbin.add.method});

    auto data = load_data(params);
    auto alignment = align(params, data, observer);
    auto geno = select(data, alignment, params.sbin.has_dom);

    notify(
        observer,
        PredictDataLoadedEvent{
            .num_samples = static_cast<size_t>(data.fam_df.rows()),
            .num_snps = static_cast<size_t>(params.snp_effects.rows()),
            .num_covar_terms = params.coefficients.names.size()});

    predict::detail::standardize_genotypes(geno, params.sbin);

    predict::SnpEffects effects{
        .add = params.add_effects, .dom = params.dom_effects};
    auto gebv = predict::detail::compute_gebv(geno, effects);
    auto covar = predict::detail::compute_covariate_effects(
        data.covariates, params.coefficients);

    auto sample_keys = data.fam_df.index().keys();
    std::vector<std::string> sample_ids(sample_keys.begin(), sample_keys.end());

    predict::PredictResult result{
        .sample_ids = std::move(sample_ids),
        .predictions = gebv.total + covar.total,
        .snp_predictions = std::move(gebv.total),
        .add_predictions = std::move(gebv.add_predictions),
        .dom_predictions = std::move(gebv.dom_predictions),
        .covar_predictions = std::move(covar.per_covariate),
        .covar_names = std::move(covar.covar_names)};

    predict::detail::PredictWriter writer(config_.output_path);
    writer.write(result);

    notify(
        observer,
        PredictResultsWrittenEvent{
            .output_path = config_.output_path.string(),
            .num_samples = result.sample_ids.size()});
}

}  // namespace gelex
