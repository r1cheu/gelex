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

#include "gelex/data/pipe/pheno.h"

#include <fmt/format.h>
#include <cstddef>
#include <optional>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/logging/pheno_event.h"
#include "gelex/infra/stats/rank_inverse_norm_transform.h"
#include "gelex/types/fixed_designs.h"

namespace gelex
{

PhenoPipe::PhenoPipe(const Config& config, PhenoObserver observer)
    : config_(config), observer_(std::move(observer))
{
    load_phenotypes();
    load_covariates();
}

auto PhenoPipe::sample_indices() const
    -> std::vector<const dataframe::Index<std::string>*>
{
    std::vector<const dataframe::Index<std::string>*> result;
    result.push_back(&phenotype_frame_->index());
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

auto PhenoPipe::load(const dataframe::Index<std::string>& common_index) -> void
{
    phenotype_frame_->gather(common_index);
    phenotype_ = phenotype_frame_->col(0).to_map<double>();

    std::optional<QuantitativeCovariate> qcov;
    std::optional<DiscreteCovariate> dcov;

    if (qcovar_frame_)
    {
        qcovar_frame_->gather(common_index);
        qcov = make_quantitative_covariate(*qcovar_frame_);
    }

    if (dcovar_frame_)
    {
        dcovar_frame_->gather(common_index);
        dcov = make_discrete_covariate(*dcovar_frame_);
    }

    if (!dcov && !qcov)
    {
        fixed_design_ = FixedDesign::make(phenotype_.size());
    }
    else
    {
        fixed_design_ = FixedDesign::make(std::move(qcov), std::move(dcov));
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

    stats::RankInverseNormTransform transformer(offset);
    auto logger = gelex::logging::get();

    if (type == detail::TransformType::DINT)
    {
        logger->info("   Method: Direct INT (DINT), offset (k): {}", offset);
        transformer.apply_dint(phenotype_);
    }
    else if (type == detail::TransformType::IINT)
    {
        logger->info("   Method: Indirect INT (IINT), offset (k): {}", offset);
        transformer.apply_iint(phenotype_, fixed_design_.X);
        fixed_design_ = FixedDesign::make(phenotype_.size());
    }
}

}  // namespace gelex
