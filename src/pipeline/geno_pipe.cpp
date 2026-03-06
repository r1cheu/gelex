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

#include <memory>
#include <utility>

#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/logging/notify.h"

namespace gelex
{

GenoPipe::GenoPipe(const Config& config, DataPipeObserver observer)
    : config_(config), observer_(std::move(observer))
{
}

auto GenoPipe::load(std::shared_ptr<SampleManager> sample_manager) -> void
{
    sample_manager_ = std::move(sample_manager);

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

auto GenoPipe::load_additive_matrix() -> void
{
    load_genotype_impl<GeneticEffectType::Add>(
        ".add", config_.genotype_method, additive_matrix_);
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

auto GenoPipe::load_dominance_matrix() -> void
{
    load_genotype_impl<GeneticEffectType::Dom>(
        ".dom", config_.genotype_method, dominance_matrix_);
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

}  // namespace gelex
