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

#include "gelex/pipeline/data_pipe.h"

#include <algorithm>
#include <filesystem>
#include <format>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/frame/dummy_encode.h"
#include "gelex/data/genotype/sample_manager.h"
#include "gelex/data/grm/grm_loader.h"
#include "gelex/exception.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/utils/formatter.h"
#include "gelex/infra/utils/phenotype_transformer.h"
#include "gelex/types/fixed_effects.h"

namespace gelex
{
DataPipe::~DataPipe() = default;
DataPipe::DataPipe(DataPipe&&) noexcept = default;
DataPipe& DataPipe::operator=(DataPipe&&) noexcept = default;

DataPipe::DataPipe(const Config& config, DataPipeObserver observer)
    : config_(config), observer_(std::move(observer))
{
    auto fam_path = config.bed_path;
    sample_manager_
        = std::make_shared<SampleManager>(fam_path.replace_extension(".fam"));
    num_genotype_samples_ = sample_manager_->num_common_samples();
}

auto DataPipe::load() -> void
{
    load_phenotypes();
    load_covariates();

    intersect_samples();

    if (!config_.grm_paths.empty())
    {
        load_grms();
    }
    else
    {
        if (config_.model_type == ModelType::A)
        {
            load_additive_matrix();
        }
        else if (config_.model_type == ModelType::D)
        {
            load_dominance_matrix();
        }
        else
        {
            load_additive_matrix();
            load_dominance_matrix();
        }
    }
    finalize();
}

auto DataPipe::load_phenotypes() -> void
{
    if (config_.phenotype_path.empty())
    {
        throw ArgumentValidationException("Phenotype file path is required.");
    }

    PhenotypeLoadedEvent event;
    event.geno_samples = num_genotype_samples_;

    int column_index
        = config_.phenotype_column - 2;  // zero-based index for data frame

    auto frame = DataFrame<double>::read(config_.phenotype_path);

    if (column_index < 0 || column_index >= static_cast<int>(frame.ncols()))
    {
        throw ColumnRangeException(
            std::format(
                "Phenotype column {} is out of range, expected [2, {}]",
                config_.phenotype_column,
                frame.ncols() + 2));
    }

    phenotype_name_
        = frame.column(static_cast<size_t>(config_.phenotype_column - 2))
              .name();

    phenotype_frame_ = std::move(frame);
    event.pheno_samples = phenotype_frame_.nrows();
    event.trait_name = phenotype_name_;

    if (observer_)
    {
        observer_(DataPipeSectionEvent{});
        observer_(event);
    }
}

auto DataPipe::load_covariates() -> void
{
    CovariatesLoadedEvent event;

    if (config_.quantitative_covar_path)
    {
        qcovar_frame_
            = DataFrame<double>::read(*config_.quantitative_covar_path);
        event.num_quantitative_covariates = qcovar_frame_->ncols();
        event.quantitative_names = qcovar_frame_->columns();
    }

    if (config_.discrete_covar_path)
    {
        dcovar_frame_
            = DataFrame<std::string>::read(*config_.discrete_covar_path);
        event.num_discrete_covariates = dcovar_frame_->ncols();
        event.discrete_names = dcovar_frame_->columns();
    }

    if (observer_
        && (event.num_quantitative_covariates || event.num_discrete_covariates))
    {
        observer_(event);
    }
}

auto DataPipe::load_grms() -> void
{
    for (const auto& grm_path : config_.grm_paths)
    {
        grm_loaders_.emplace_back(grm_path);
        if (observer_)
        {
            observer_(
                GrmLoadedEvent{
                    .num_samples
                    = static_cast<size_t>(grm_loaders_.back().num_samples()),
                    .type = grm_loaders_.back().type()});
        }
    }
}

auto DataPipe::intersect_samples() -> void
{
    size_t total_before = sample_manager_->num_common_samples();

    if (phenotype_frame_.nrows() == 0)
    {
        throw InvalidOperationException(
            "Phenotype frame cannot be empty."
            " Load a non-empty phenotype file first.");
    }

    total_before = std::max(total_before, phenotype_frame_.nrows());
    sample_manager_->intersect(phenotype_frame_.index_column().data());

    if (qcovar_frame_)
    {
        total_before = std::max(total_before, qcovar_frame_->nrows());
        sample_manager_->intersect(qcovar_frame_->index_column().data());
    }

    if (dcovar_frame_)
    {
        total_before = std::max(total_before, dcovar_frame_->nrows());
        sample_manager_->intersect(dcovar_frame_->index_column().data());
    }

    for (const auto& grm_loader : grm_loaders_)
    {
        total_before = std::max(
            total_before, static_cast<size_t>(grm_loader.num_samples()));
        sample_manager_->intersect(grm_loader.sample_ids());
    }

    sample_manager_->finalize();

    size_t common = sample_manager_->num_common_samples();

    if (common == 0)
    {
        throw InvalidInputException(
            "No common samples found between phenotype, covariates, and "
            "genotype");
    }

    if (observer_)
    {
        observer_(
            IntersectionEvent{
                .common_samples = common,
                .excluded_samples = total_before - common});
    }
}

auto DataPipe::load_additive_matrix() -> void
{
    load_genotype_impl<GeneticEffectType::Add>(
        ".add", config_.genotype_method, additive_matrix_);
    if (observer_)
    {
        int64_t mono = 0;
        int64_t total = 0;
        std::visit(
            [&](auto&& m)
            {
                mono = m.num_mono();
                total = m.cols();
            },
            *additive_matrix_);
        observer_(
            GenotypeLoadedEvent{
                .is_dominance = false,
                .num_snps = total,
                .monomorphic_snps = mono});
    }
}

auto DataPipe::load_dominance_matrix() -> void
{
    load_genotype_impl<GeneticEffectType::Dom>(
        ".dom", config_.genotype_method, dominance_matrix_);
    if (observer_)
    {
        int64_t mono = 0;
        int64_t total = 0;
        std::visit(
            [&](auto&& m)
            {
                mono = m.num_mono();
                total = m.cols();
            },
            *dominance_matrix_);
        observer_(
            GenotypeLoadedEvent{
                .is_dominance = true,
                .num_snps = total,
                .monomorphic_snps = mono});
    }
}

auto DataPipe::finalize() -> void
{
    const auto& common_ids = sample_manager_->common_ids();
    const auto& id_map = sample_manager_->common_id_map();

    if (phenotype_frame_.nrows() == 0)
    {
        throw InvalidOperationException(
            "Phenotype frame cannot be empty."
            " Load a non-empty phenotype file first.");
    }

    auto aligned = phenotype_frame_;
    aligned.intersect_index_inplace(common_ids);

    const auto& values
        = aligned.column(static_cast<size_t>(config_.phenotype_column - 2))
              .data();
    phenotype_ = Eigen::Map<const Eigen::VectorXd>(
        values.data(), static_cast<Eigen::Index>(values.size()));

    std::optional<QuantitativeCovariate> qcov;
    std::optional<DiscreteCovariate> dcov;

    if (qcovar_frame_)
    {
        auto aligned = *qcovar_frame_;
        aligned.intersect_index_inplace(common_ids);

        auto qcov_matrix = aligned.eigen();
        auto names = aligned.columns();

        qcov = QuantitativeCovariate{
            .names = std::move(names), .X = std::move(qcov_matrix)};
    }

    if (dcovar_frame_)
    {
        auto aligned = *dcovar_frame_;
        aligned.intersect_index_inplace(common_ids);

        dcov = DummyEncode(aligned);
    }
    if (!dcov && !qcov)
    {
        fixed_effects_ = FixedEffect::build(phenotype_.size());
    }
    else
    {
        fixed_effects_ = FixedEffect::build(std::move(qcov), std::move(dcov));
    }
    grms_.reserve(grm_loaders_.size());
    for (auto& grm_loader : grm_loaders_)
    {
        grms_.emplace_back(grm_loader.type(), grm_loader.load(id_map));
    }

    apply_phenotype_transform(config_.transform_type, config_.int_offset);
}

auto DataPipe::apply_phenotype_transform(
    detail::TransformType type,
    double offset) -> void
{
    if (type == detail::TransformType::None)
    {
        return;
    }

    detail::PhenotypeTransformer transformer(offset);
    auto logger = gelex::logging::get();

    if (type == detail::TransformType::DINT)
    {
        logger->info(task("Method: Direct INT (DINT), offset (k): {}", offset));
        transformer.apply_dint(phenotype_);
    }
    else if (type == detail::TransformType::IINT)
    {
        logger->info(
            task("Method: Indirect INT (IINT), offset (k): {}", offset));
        transformer.apply_iint(phenotype_, fixed_effects_.X);
        fixed_effects_ = FixedEffect::build(phenotype_.size());
    }
}

}  // namespace gelex
