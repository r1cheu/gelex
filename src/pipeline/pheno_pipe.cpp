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
#include <algorithm>
#include <optional>
#include <span>
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

auto PhenoPipe::gather_by_ids(
    df::DataFrame<std::string>& frame,
    std::span<const std::string> ids) -> void
{
    std::vector<std::size_t> pos;
    pos.reserve(ids.size());
    for (const auto& id : ids)
    {
        pos.push_back(frame.index().at(id));
    }
    frame.gather(pos);
}

auto PhenoPipe::build_discrete_covariate(
    const df::DataFrame<std::string>& frame) -> DiscreteCovariate
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

PhenoPipe::PhenoPipe(const Config& config, DataPipeObserver observer)
    : config_(config), observer_(std::move(observer))
{
    auto fam_path = config.bed_path;
    auto fam = read_fam(fam_path.replace_extension(".fam"));
    sample_index_ = std::move(fam).index();
    num_genotype_samples_ = sample_index_.size();

    notify(observer_, DataPipeSectionEvent{});
}

auto PhenoPipe::load(const std::vector<std::span<const std::string>>& extra_ids)
    -> void
{
    load_phenotypes();
    load_covariates();
    intersect_samples(extra_ids);
    finalize();
}

auto PhenoPipe::load_phenotypes() -> void
{
    if (config_.phenotype_path.empty())
    {
        throw GelexException("Phenotype file path is required.");
    }

    PhenotypeLoadedEvent event;
    event.geno_samples = num_genotype_samples_;

    if (config_.phenotype_column < 2)
    {
        throw GelexException(
            fmt::format(
                "Phenotype column {} is out of range, expected >= 2",
                config_.phenotype_column));
    }

    auto pheno_col = static_cast<std::size_t>(config_.phenotype_column - 2);
    auto frame = read_pheno(config_.phenotype_path, &pheno_col);

    phenotype_name_ = std::string(frame.col(0).name());

    phenotype_frame_ = std::move(frame);
    event.pheno_samples = phenotype_frame_->rows();
    event.trait_name = phenotype_name_;

    notify(observer_, event);
}

auto PhenoPipe::load_covariates() -> void
{
    CovariatesLoadedEvent event;

    if (config_.quantitative_covariates_path)
    {
        qcovar_frame_ = read_qcovar(*config_.quantitative_covariates_path);
        event.num_quantitative_covariates = qcovar_frame_->cols();
        event.quantitative_names = std::vector<std::string>(
            qcovar_frame_->names().begin(), qcovar_frame_->names().end());
    }

    if (config_.discrete_covariates_path)
    {
        dcovar_frame_ = read_dcovar(*config_.discrete_covariates_path);
        event.num_discrete_covariates = dcovar_frame_->cols();
        event.discrete_names = std::vector<std::string>(
            dcovar_frame_->names().begin(), dcovar_frame_->names().end());
    }

    if (event.num_quantitative_covariates || event.num_discrete_covariates)
    {
        notify(observer_, event);
    }
}

auto PhenoPipe::intersect_samples(
    const std::vector<std::span<const std::string>>& extra_ids) -> void
{
    if (!phenotype_frame_ || phenotype_frame_->rows() == 0)
    {
        throw GelexException(
            "Phenotype frame cannot be empty."
            " Load a non-empty phenotype file first.");
    }

    std::vector<df::Index<std::string>> extra_indices;
    extra_indices.reserve(extra_ids.size());
    for (const auto& ids : extra_ids)
    {
        extra_indices.emplace_back(
            std::vector<std::string>(ids.begin(), ids.end()));
    }

    std::vector<const df::Index<std::string>*> all_indices;
    all_indices.push_back(&sample_index_);
    all_indices.push_back(&phenotype_frame_->index());
    if (qcovar_frame_)
    {
        all_indices.push_back(&qcovar_frame_->index());
    }
    if (dcovar_frame_)
    {
        all_indices.push_back(&dcovar_frame_->index());
    }
    for (const auto& idx : extra_indices)
    {
        all_indices.push_back(&idx);
    }

    size_t total_before = 0;
    for (const auto* idx : all_indices)
    {
        total_before = std::max(total_before, idx->size());
    }

    auto positions = df::intersect<std::string>(all_indices);
    sample_index_.gather(positions[0]);

    size_t common = sample_index_.size();

    if (common == 0)
    {
        throw GelexException(
            "No common samples found between phenotype, covariates, and "
            "genotype");
    }

    notify(
        observer_,
        IntersectionEvent{
            .common_samples = common,
            .excluded_samples = total_before - common});
}

auto PhenoPipe::finalize() -> void
{
    const auto common_ids = sample_index_.keys();

    if (!phenotype_frame_ || phenotype_frame_->rows() == 0)
    {
        throw GelexException(
            "Phenotype frame cannot be empty."
            " Load a non-empty phenotype file first.");
    }

    auto aligned = phenotype_frame_->clone();
    gather_by_ids(aligned, common_ids);

    auto values = aligned.col(0).template as<double>();
    phenotype_ = Eigen::Map<const Eigen::VectorXd>(
        values.data(), static_cast<Eigen::Index>(values.size()));

    std::optional<QuantitativeCovariate> qcov;
    std::optional<DiscreteCovariate> dcov;

    if (qcovar_frame_)
    {
        auto qcovar_aligned = qcovar_frame_->clone();
        gather_by_ids(qcovar_aligned, common_ids);

        auto names = std::vector<std::string>(
            qcovar_aligned.names().begin(), qcovar_aligned.names().end());
        qcov = QuantitativeCovariate{
            .names = std::move(names), .X = qcovar_aligned.to_mat<double>()};
    }

    if (dcovar_frame_)
    {
        auto dcovar_aligned = dcovar_frame_->clone();
        gather_by_ids(dcovar_aligned, common_ids);
        dcov = build_discrete_covariate(dcovar_aligned);
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
