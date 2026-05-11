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

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype/genotype.h"
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
        std::string bfile_prefix;

        std::vector<GeneticMode> requested_effects;
        GenotypeProcessMethod genotype_method;
        bool use_mmap = false;
        int chunk_size = 10000;

        std::string output_prefix;
    };

    explicit GenoPipe(const Config& config, GenoObserver observer = {});
    GenoPipe(const GenoPipe&) = delete;
    GenoPipe(GenoPipe&&) noexcept = default;
    auto operator=(const GenoPipe&) -> GenoPipe& = delete;
    auto operator=(GenoPipe&&) noexcept -> GenoPipe& = default;
    ~GenoPipe() = default;

    auto load(const dataframe::Index<std::string>& sample_index) -> void;

    auto take_additive_matrix() && -> genotype::Genotype
    {
        return std::move(*additive_matrix_);
    }

    auto take_dominance_matrix() && -> genotype::Genotype
    {
        return std::move(*dominance_matrix_);
    }

    [[nodiscard]] auto has_additive_matrix() const -> bool
    {
        return additive_matrix_.has_value();
    }

    [[nodiscard]] auto has_dominance_matrix() const -> bool
    {
        return dominance_matrix_.has_value();
    }

    [[nodiscard]] auto sample_indices() const
        -> std::vector<const dataframe::Index<std::string>*>
    {
        return {&fam_index_};
    }

   private:
    template <GeneticMode GT>
    auto load_genotype_impl(
        const dataframe::Index<std::string>& sample_index,
        const std::string& suffix,
        GenotypeProcessMethod method,
        std::optional<genotype::Genotype>& target) -> void;

    auto load_additive_matrix(const dataframe::Index<std::string>& sample_index)
        -> void;
    auto load_dominance_matrix(
        const dataframe::Index<std::string>& sample_index) -> void;
    auto write_sbin() -> void;

    Config config_;
    dataframe::Index<std::string> fam_index_;
    std::optional<genotype::Genotype> additive_matrix_;
    std::optional<genotype::Genotype> dominance_matrix_;
    GenoObserver observer_;
};

}  // namespace gelex

#endif  // GELEX_DATA_PIPE_GENO_H_
