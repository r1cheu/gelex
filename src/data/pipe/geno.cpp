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

#include "gelex/data/pipe/geno.h"

#include <Eigen/Core>
#include <algorithm>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype/genotype.h"
#include "gelex/data/genotype/genotype_reader.h"
#include "gelex/data/genotype/method.h"
#include "gelex/data/reader.h"
#include "gelex/infra/logging/geno_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/locistats/writer.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

GenoPipe::GenoPipe(const Config& config, GenoObserver observer)
    : config_(config),
      fam_index_(read_fam(config_.bfile_prefix + ".fam").index()),
      observer_(std::move(observer))
{
}

auto GenoPipe::load(const dataframe::Index<std::string>& sample_index) -> void
{
    if (std::ranges::contains(config_.requested_effects, GeneticMode::A))
    {
        load_additive_matrix(sample_index);
    }
    if (std::ranges::contains(config_.requested_effects, GeneticMode::D))
    {
        load_dominance_matrix(sample_index);
    }
    write_sbin();
}

template <GeneticMode GT>
auto GenoPipe::load_genotype_impl(
    const dataframe::Index<std::string>& sample_index,
    const std::string& suffix,
    GenotypeMethod method,
    std::optional<genotype::Genotype>& target) -> void
{
    genotype::GenotypeReader reader(
        config_.bfile_prefix, sample_index, observer_);
    if (config_.use_mmap)
    {
        std::filesystem::path prefix = config_.output_prefix + suffix;
        target.emplace(reader.template read<GT>(
            method,
            genotype::GenotypeReader::Sink::Mmap{std::move(prefix)},
            config_.chunk_size));
    }
    else
    {
        target.emplace(reader.template read<GT>(
            method,
            genotype::GenotypeReader::Sink::InMemory{},
            config_.chunk_size));
    }
}

auto GenoPipe::load_additive_matrix(
    const dataframe::Index<std::string>& sample_index) -> void
{
    load_genotype_impl<GeneticMode::A>(
        sample_index, ".add", config_.genotype_method, additive_matrix_);
    notify(
        observer_,
        GenotypeLoadedEvent{
            .mode = GeneticMode::A,
            .num_snps = additive_matrix_->cols(),
            .monomorphic_snps = additive_matrix_->num_mono()});
}

auto GenoPipe::load_dominance_matrix(
    const dataframe::Index<std::string>& sample_index) -> void
{
    load_genotype_impl<GeneticMode::D>(
        sample_index, ".dom", config_.genotype_method, dominance_matrix_);
    notify(
        observer_,
        GenotypeLoadedEvent{
            .mode = GeneticMode::D,
            .num_snps = dominance_matrix_->cols(),
            .monomorphic_snps = dominance_matrix_->num_mono()});
}

auto GenoPipe::write_sbin() -> void
{
    LociStatsWriter writer(config_.output_prefix + ".sbin");
    auto method_code = std::to_underlying(config_.genotype_method);
    const bool method_is_center = is_center(config_.genotype_method);

    if (additive_matrix_)
    {
        writer.write(
            EffectType::add(),
            method_code,
            additive_matrix_->mean(),
            method_is_center ? static_cast<const Eigen::VectorXd*>(nullptr)
                             : &additive_matrix_->stddev(),
            additive_matrix_->mono_indices());
    }

    if (dominance_matrix_)
    {
        writer.write(
            EffectType::dom(),
            method_code,
            dominance_matrix_->mean(),
            method_is_center ? static_cast<const Eigen::VectorXd*>(nullptr)
                             : &dominance_matrix_->stddev(),
            dominance_matrix_->mono_indices());
    }
}

}  // namespace gelex
