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

#include "gelex/pipeline/geno_pipe.h"

#include <utility>

#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/locistats_writer.h"
#include "gelex/types/genotype_process_method.h"

namespace gelex
{

GenoPipe::GenoPipe(const Config& config, DataPipeObserver observer)
    : config_(config), observer_(std::move(observer))
{
}

auto GenoPipe::load(const df::Index<std::string>& sample_index) -> void
{
    if (config_.model_type == ModelType::A)
    {
        load_additive_matrix(sample_index);
    }
    else if (config_.model_type == ModelType::D)
    {
        load_dominance_matrix(sample_index);
    }
    else
    {
        load_additive_matrix(sample_index);
        load_dominance_matrix(sample_index);
    }

    write_sbin();
}

auto GenoPipe::load_additive_matrix(const df::Index<std::string>& sample_index)
    -> void
{
    load_genotype_impl<GeneticKind::Add>(
        sample_index, ".add", config_.genotype_method, additive_matrix_);
    int64_t mono = 0;
    int64_t total = 0;
    std::visit(
        [&](auto&& m)
        {
            mono = m.num_mono();
            total = m.cols();
        },
        *additive_matrix_);
    notify(
        observer_,
        GenotypeLoadedEvent{
            .is_dominance = false,
            .num_snps = total,
            .monomorphic_snps = mono});
}

auto GenoPipe::load_dominance_matrix(const df::Index<std::string>& sample_index)
    -> void
{
    load_genotype_impl<GeneticKind::Dom>(
        sample_index, ".dom", config_.genotype_method, dominance_matrix_);
    int64_t mono = 0;
    int64_t total = 0;
    std::visit(
        [&](auto&& m)
        {
            mono = m.num_mono();
            total = m.cols();
        },
        *dominance_matrix_);
    notify(
        observer_,
        GenotypeLoadedEvent{
            .is_dominance = true, .num_snps = total, .monomorphic_snps = mono});
}

auto GenoPipe::write_sbin() -> void
{
    LociStatsWriter writer(config_.output_prefix + ".sbin");
    auto method_code = static_cast<uint8_t>(config_.genotype_method);
    const bool is_center = is_center_family_method(config_.genotype_method);

    if (additive_matrix_)
    {
        std::visit(
            [&](const auto& m)
            {
                writer.write(
                    EffectType::add(),
                    method_code,
                    m.mean(),
                    is_center ? static_cast<const Eigen::VectorXd*>(nullptr)
                              : &m.stddev(),
                    m.mono_indices());
            },
            *additive_matrix_);
    }

    if (dominance_matrix_)
    {
        std::visit(
            [&](const auto& m)
            {
                writer.write(
                    EffectType::dom(),
                    method_code,
                    m.mean(),
                    is_center ? static_cast<const Eigen::VectorXd*>(nullptr)
                              : &m.stddev(),
                    m.mono_indices());
            },
            *dominance_matrix_);
    }

    writer.finalize();
}

}  // namespace gelex
