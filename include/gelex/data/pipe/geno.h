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

#ifndef GELEX_DATA_PIPE_GENO_H_
#define GELEX_DATA_PIPE_GENO_H_

#include <filesystem>
#include <string>
#include <variant>
#include <vector>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype/map_reader.h"
#include "gelex/data/genotype/mat_reader.h"
#include "gelex/data/genotype/matrix.h"
#include "gelex/data/genotype/mmap.h"
#include "gelex/data/genotype/process_method.h"
#include "gelex/infra/logging/geno_event.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class GenoPipe
{
   public:
    struct Config
    {
        std::filesystem::path bed_path;

        GeneticMode mode;
        GenotypeProcessMethod genotype_method;
        bool use_mmap = false;
        int chunk_size = 10000;

        std::string output_prefix;
    };

    explicit GenoPipe(const Config& config, GenoObserver observer = {});
    GenoPipe(const GenoPipe&) = delete;
    GenoPipe(GenoPipe&&) noexcept = default;
    GenoPipe& operator=(const GenoPipe&) = delete;
    GenoPipe& operator=(GenoPipe&&) noexcept = default;
    ~GenoPipe() = default;

    auto load(const dataframe::Index<std::string>& sample_index) -> void;

    auto take_additive_matrix() && -> std::
        variant<genotype::GenotypeMap, genotype::GenotypeMatrix>
    {
        return std::move(*additive_matrix_);
    }

    auto take_dominance_matrix() && -> std::
        variant<genotype::GenotypeMap, genotype::GenotypeMatrix>
    {
        return std::move(*dominance_matrix_);
    }

    auto has_additive_matrix() const -> bool
    {
        return additive_matrix_ != nullptr;
    }

    auto has_dominance_matrix() const -> bool
    {
        return dominance_matrix_ != nullptr;
    }

    auto sample_indices() const
        -> std::vector<const dataframe::Index<std::string>*>
    {
        return {&fam_index_};
    }

   private:
    using GenotypeMatrixPtr = std::unique_ptr<
        std::variant<genotype::GenotypeMap, genotype::GenotypeMatrix>>;

    template <GeneticMode GT>
    auto load_genotype_impl(
        const dataframe::Index<std::string>& sample_index,
        const std::string& suffix,
        GenotypeProcessMethod method,
        GenotypeMatrixPtr& target) -> void
    {
        if (config_.use_mmap)
        {
            std::string file_path = config_.output_prefix + suffix;
            auto pipe = genotype::GenotypeMapReader(
                config_.bed_path, sample_index, file_path, observer_);
            target = std::make_unique<
                std::variant<genotype::GenotypeMap, genotype::GenotypeMatrix>>(
                pipe.process<GT>(method, config_.chunk_size));
        }
        else
        {
            auto reader = genotype::GenotypeMatReader(
                config_.bed_path, sample_index, observer_);
            target = std::make_unique<
                std::variant<genotype::GenotypeMap, genotype::GenotypeMatrix>>(
                reader.process<GT>(method, config_.chunk_size));
        }
    }

    auto load_additive_matrix(const dataframe::Index<std::string>& sample_index)
        -> void;
    auto load_dominance_matrix(
        const dataframe::Index<std::string>& sample_index) -> void;
    auto write_sbin() -> void;

    Config config_;
    dataframe::Index<std::string> fam_index_;
    GenotypeMatrixPtr additive_matrix_;
    GenotypeMatrixPtr dominance_matrix_;
    GenoObserver observer_;
};

}  // namespace gelex

#endif  // GELEX_DATA_PIPE_GENO_H_
