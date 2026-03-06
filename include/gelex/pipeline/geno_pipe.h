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

#ifndef GELEX_PIPELINE_GENO_PIPE_H_
#define GELEX_PIPELINE_GENO_PIPE_H_

#include <filesystem>
#include <memory>
#include <string>
#include <variant>

#include "gelex/data/genotype/genotype_loader.h"
#include "gelex/data/genotype/genotype_matrix.h"
#include "gelex/data/genotype/genotype_mmap.h"
#include "gelex/data/genotype/genotype_pipe.h"
#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class SampleManager;

class GenoPipe
{
   public:
    struct Config
    {
        std::filesystem::path bed_path;

        ModelType model_type;
        GenotypeProcessMethod genotype_method;
        bool use_mmap = false;
        int chunk_size = 10000;

        std::string output_prefix;
    };

    explicit GenoPipe(const Config& config, DataPipeObserver observer = {});
    GenoPipe(const GenoPipe&) = delete;
    GenoPipe(GenoPipe&&) noexcept = default;
    GenoPipe& operator=(const GenoPipe&) = delete;
    GenoPipe& operator=(GenoPipe&&) noexcept = default;
    ~GenoPipe() = default;

    auto load(std::shared_ptr<SampleManager> sample_manager) -> void;

    auto take_additive_matrix() && -> std::variant<GenotypeMap, GenotypeMatrix>
    {
        return std::move(*additive_matrix_);
    }

    auto take_dominance_matrix() && -> std::variant<GenotypeMap, GenotypeMatrix>
    {
        return std::move(*dominance_matrix_);
    }

    auto has_dominance_matrix() const -> bool
    {
        return dominance_matrix_ != nullptr;
    }

   private:
    using GenotypeMatrixPtr
        = std::unique_ptr<std::variant<GenotypeMap, GenotypeMatrix>>;

    template <GeneticEffectType GT>
    auto load_genotype_impl(
        const std::string& suffix,
        GenotypeProcessMethod method,
        GenotypeMatrixPtr& target) -> void
    {
        if (config_.use_mmap)
        {
            std::string file_path = config_.output_prefix + suffix;
            auto pipe = gelex::GenotypePipe(
                config_.bed_path, sample_manager_, file_path);
            target
                = std::make_unique<std::variant<GenotypeMap, GenotypeMatrix>>(
                    pipe.process<GT>(method, config_.chunk_size));
        }
        else
        {
            auto loader
                = gelex::GenotypeLoader(config_.bed_path, sample_manager_);
            target
                = std::make_unique<std::variant<GenotypeMap, GenotypeMatrix>>(
                    loader.process<GT>(method, config_.chunk_size));
        }
    }

    auto load_additive_matrix() -> void;
    auto load_dominance_matrix() -> void;

    Config config_;
    std::shared_ptr<SampleManager> sample_manager_;
    GenotypeMatrixPtr additive_matrix_;
    GenotypeMatrixPtr dominance_matrix_;
    DataPipeObserver observer_;
};

}  // namespace gelex

#endif  // GELEX_PIPELINE_GENO_PIPE_H_
