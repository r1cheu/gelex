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

#include "gelex/pipeline/pheno_pipe.h"

#include <fmt/format.h>
#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/dataframe/encode.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/stats/rank_inverse_norm_transform.h"
#include "gelex/types/fixed_effects.h"

namespace gelex
{

namespace
{

auto build_discrete_covariate(const df::DataFrame<std::string>& frame)
    -> DiscreteCovariate
{
    std::vector<std::string> names;
    std::vector<std::vector<std::string>> levels;
    std::vector<std::string> reference_levels;
    std::vector<df::EncodedResult<>> encoded_results;

    for (std::size_t i = 0; i < frame.cols(); ++i)
    {
        const auto& col = frame.col(i);
        auto all_levels = df::collect_levels(col);
        if (all_levels.size() < 2)
        {
            continue;
        }
        names.emplace_back(col.name());
        reference_levels.push_back(all_levels.front());
        levels.push_back(all_levels);
        encoded_results.push_back(df::dummy_encode(col));
    }

    Eigen::Index total_cols = 0;
    for (const auto& r : encoded_results)
    {
        total_cols += r.data.cols();
    }

    Eigen::MatrixXd X(static_cast<Eigen::Index>(frame.rows()), total_cols);
    Eigen::Index col_offset = 0;
    for (const auto& r : encoded_results)
    {
        X.middleCols(col_offset, r.data.cols()) = r.data;
        col_offset += r.data.cols();
    }

    return DiscreteCovariate{
        .names = std::move(names),
        .levels = std::move(levels),
        .reference_levels = std::move(reference_levels),
        .X = std::move(X)};
}

}  // namespace

PhenoPipe::PhenoPipe(const Config& config, DataPipeObserver observer)
    : config_(config), observer_(std::move(observer))
{
    notify(observer_, DataPipeSectionEvent{});
}

auto PhenoPipe::load() -> void
{
    load_phenotypes();
    load_covariates();
}

auto PhenoPipe::covar_indices() const
    -> std::vector<const df::Index<std::string>*>
{
    std::vector<const df::Index<std::string>*> result;
    if (qcovar_frame_)
    {
        result.push_back(&qcovar_frame_->index());
    }
    if (dcovar_frame_)
    {
        result.push_back(&dcovar_frame_->index());
    }
    return result;
}

auto PhenoPipe::load_phenotypes() -> void
{
    if (config_.phenotype_path.empty())
    {
        throw GelexException("Phenotype file path is required.");
    }

    if (config_.phenotype_column < 2)
    {
        throw GelexException(
            fmt::format(
                "Phenotype column {} is out of range, expected >= 2",
                config_.phenotype_column));
    }

    auto pheno_col = static_cast<std::size_t>(config_.phenotype_column - 2);
    auto frame = read_pheno(config_.phenotype_path, &pheno_col);

    auto trait_name = std::string(frame.col(0).name());

    phenotype_frame_ = std::move(frame);

    notify(
        observer_,
        PhenotypeLoadedEvent{
            .pheno_samples = phenotype_frame_->rows(),
            .trait_name = std::move(trait_name)});
}

auto PhenoPipe::load_covariates() -> void
{
    CovariatesLoadedEvent event;

    if (config_.quantitative_covariates_path)
    {
        qcovar_frame_ = read_qcovar(*config_.quantitative_covariates_path);
        event.num_quantitative_covariates = qcovar_frame_->cols();
        event.quantitative_names
            = std::ranges::to<std::vector<std::string>>(qcovar_frame_->names());
    }

    if (config_.discrete_covariates_path)
    {
        dcovar_frame_ = read_dcovar(*config_.discrete_covariates_path);
        event.num_discrete_covariates = dcovar_frame_->cols();
        event.discrete_names
            = std::ranges::to<std::vector<std::string>>(dcovar_frame_->names());
    }

    if (event.num_quantitative_covariates || event.num_discrete_covariates)
    {
        notify(observer_, event);
    }
}

auto PhenoPipe::gather(const df::Index<std::string>& common_index) -> void
{
    phenotype_frame_->gather(common_index);
    phenotype_ = phenotype_frame_->col(0).to_map<double>();

    std::optional<QuantitativeCovariate> qcov;
    std::optional<DiscreteCovariate> dcov;

    if (qcovar_frame_)
    {
        qcovar_frame_->gather(common_index);
        auto names = qcovar_frame_->names();
        qcov = QuantitativeCovariate{
            .names = std::ranges::to<std::vector<std::string>>(names),
            .X = qcovar_frame_->to_mat<double>()};
    }

    if (dcovar_frame_)
    {
        dcovar_frame_->gather(common_index);
        dcov = build_discrete_covariate(*dcovar_frame_);
    }

    if (!dcov && !qcov)
    {
        fixed_effects_ = FixedEffect::build(phenotype_.size());
    }
    else
    {
        fixed_effects_ = FixedEffect::build(std::move(qcov), std::move(dcov));
    }

    apply_phenotype_transform(config_.transform_type, config_.int_offset);
}

auto PhenoPipe::apply_phenotype_transform(
    detail::TransformType type,
    double offset) -> void
{
    if (type == detail::TransformType::None)
    {
        return;
    }

    RankInverseNormTransform transformer(offset);
    auto logger = gelex::logging::get();

    if (type == detail::TransformType::DINT)
    {
        logger->info("   Method: Direct INT (DINT), offset (k): {}", offset);
        transformer.apply_dint(phenotype_);
    }
    else if (type == detail::TransformType::IINT)
    {
        logger->info("   Method: Indirect INT (IINT), offset (k): {}", offset);
        transformer.apply_iint(phenotype_, fixed_effects_.X);
        fixed_effects_ = FixedEffect::build(phenotype_.size());
    }
}

}  // namespace gelex
