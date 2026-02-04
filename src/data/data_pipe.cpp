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

#include "gelex/data/data_pipe.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "../src/utils/phenotype_transformer.h"
#include "gelex/data/genotype_processor.h"
#include "gelex/data/sample_manager.h"
#include "gelex/exception.h"
#include "gelex/logger.h"
#include "grm_loader.h"
#include "loader/dcovariate_loader.h"
#include "loader/phenotype_loader.h"
#include "loader/qcovariate_loader.h"
#include "types/fixed_effects.h"
#include "utils/formatter.h"

namespace bk = barkeep;

namespace gelex
{
DataPipe::~DataPipe() = default;
DataPipe::DataPipe(DataPipe&&) noexcept = default;
DataPipe& DataPipe::operator=(DataPipe&&) noexcept = default;

DataPipe::DataPipe(const Config& config) : config_(config)
{
    auto fam_path = config.bed_path;
    fam_path.replace_extension(".fam");
    sample_manager_ = std::make_shared<SampleManager>(fam_path);
    num_genotype_samples_ = sample_manager_->num_common_samples();
}

PhenoStats DataPipe::load_phenotypes()
{
    if (config_.phenotype_path.empty())
    {
        throw ArgumentValidationException("Phenotype file path is required.");
    }

    phenotype_loader_ = std::make_unique<detail::PhenotypeLoader>(
        config_.phenotype_path, config_.phenotype_column, config_.iid_only);

    return PhenoStats{
        .samples_loaded = phenotype_loader_->sample_ids().size(),
        .trait_name = phenotype_loader_->name()};
}

CovarStats DataPipe::load_covariates()
{
    std::vector<std::string> q_names;
    std::vector<std::string> d_names;

    if (!config_.qcovar_path.empty())
    {
        qcovar_loader_ = std::make_unique<detail::QuantitativeCovariateLoader>(
            config_.qcovar_path, config_.iid_only);
        q_names = qcovar_loader_->column_names();
    }

    if (!config_.dcovar_path.empty())
    {
        dcovar_loader_ = std::make_unique<detail::DiscreteCovariateLoader>(
            config_.dcovar_path, config_.iid_only);
        d_names = dcovar_loader_->column_names();
    }

    return CovarStats{
        .qcovar_loaded = q_names.size(),
        .dcovar_loaded = d_names.size(),
        .q_names = std::move(q_names),
        .d_names = std::move(d_names)};
}

auto DataPipe::load_grms() -> std::vector<GrmStats>
{
    std::vector<GrmStats> grm_stats;
    grm_stats.reserve(config_.grm_paths.size());

    for (const auto& grm_path : config_.grm_paths)
    {
        grm_loaders_.emplace_back(grm_path);
        grm_stats.push_back(
            GrmStats{
                .samples_in_file
                = static_cast<size_t>(grm_loaders_.back().num_samples()),
                .type = grm_loaders_.back().type()});
    }
    return grm_stats;
}

IntersectionStats DataPipe::intersect_samples()
{
    size_t total_before = sample_manager_->num_common_samples();

    if (phenotype_loader_)
    {
        total_before
            = std::max(total_before, phenotype_loader_->sample_ids().size());
        sample_manager_->intersect(phenotype_loader_->sample_ids());
    }

    if (qcovar_loader_)
    {
        total_before
            = std::max(total_before, qcovar_loader_->sample_ids().size());
        sample_manager_->intersect(qcovar_loader_->sample_ids());
    }

    if (dcovar_loader_)
    {
        total_before
            = std::max(total_before, dcovar_loader_->sample_ids().size());
        sample_manager_->intersect(dcovar_loader_->sample_ids());
    }

    for (const auto& grm_loader : grm_loaders_)
    {
        total_before = std::max(
            total_before, static_cast<size_t>(grm_loader.num_samples()));
        sample_manager_->intersect(grm_loader.sample_ids());
    }

    sample_manager_->finalize();

    size_t common = sample_manager_->num_common_samples();

    return IntersectionStats{
        .total_samples = total_before,
        .common_samples = common,
        .excluded_samples = total_before - common};
}

GenotypeStats DataPipe::load_additive_matrix()
{
    load_genotype_impl<gelex::AdditiveProcessor<OrthStandardizeMethod>>(
        ".add", additive_matrix_);

    int64_t num_mono_snps = 0;
    int64_t total_snps = 0;

    if (additive_matrix_)
    {
        std::visit(
            [&](auto&& arg)
            {
                num_mono_snps = arg.num_mono();
                total_snps = arg.cols();
            },
            *additive_matrix_);
    }

    return GenotypeStats{
        .num_snps = total_snps, .monomorphic_snps = num_mono_snps};
}

GenotypeStats DataPipe::load_dominance_matrix()
{
    load_genotype_impl<gelex::DominantProcessor<OrthStandardizeMethod>>(
        ".dom", dominance_matrix_);

    int64_t num_mono_snps = 0;
    int64_t total_snps = 0;

    if (dominance_matrix_)
    {
        std::visit(
            [&](auto&& arg)
            {
                num_mono_snps = arg.num_mono();
                total_snps = arg.cols();
            },
            *dominance_matrix_);
    }

    return GenotypeStats{
        .num_snps = total_snps, .monomorphic_snps = num_mono_snps};
}

void DataPipe::finalize()
{
    const auto& id_map = sample_manager_->common_id_map();

    phenotype_ = phenotype_loader_ ? phenotype_loader_->load(id_map)
                                   : Eigen::VectorXd::Zero();

    std::optional<QuantitativeCovariate> qcov;
    std::optional<DiscreteCovariate> dcov;

    if (qcovar_loader_)
    {
        qcov = qcovar_loader_->load(id_map);
    }

    if (dcovar_loader_)
    {
        dcov = dcovar_loader_->load(id_map);
    }
    if (!dcov && !qcov)
    {
        fixed_effects_ = FixedEffect::build_bayes(phenotype_.size());
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

const std::vector<std::string>& DataPipe::fixed_effect_names() const
{
    return fixed_effects_.names;
}

void DataPipe::apply_phenotype_transform(
    detail::TransformType type,
    double offset)
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
