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

#ifndef GELEX_DATA_GENOTYPE_READER_H_
#define GELEX_DATA_GENOTYPE_READER_H_

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>

#include <Eigen/Core>

#include "gelex/data/bed.h"
#include "gelex/data/genotype.h"
#include "gelex/data/genotype_method.h"
#include "gelex/infra/logging/geno_event.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::genotype
{

class GenotypeReader
{
   public:
    GenotypeReader(
        const std::string& bfile_prefix,
        const dataframe::Index<std::string>& sample_index,
        gelex::GenoObserver observer = {});

    GenotypeReader(const GenotypeReader&) = delete;
    auto operator=(const GenotypeReader&) -> GenotypeReader& = delete;
    GenotypeReader(GenotypeReader&&) noexcept = default;
    auto operator=(GenotypeReader&&) noexcept -> GenotypeReader& = default;
    ~GenotypeReader() = default;

    auto read_in_memory(
        gelex::GeneticMode mode,
        gelex::GenotypeMethod method,
        std::size_t chunk_size = 10000) -> Genotype;

    auto read_mmap(
        gelex::GeneticMode mode,
        gelex::GenotypeMethod method,
        const std::filesystem::path& output_prefix,
        std::size_t chunk_size = 10000) -> Genotype;

    static auto read(
        const std::filesystem::path& gbin_path,
        gelex::GeneticMode mode) -> Genotype;

   private:
    struct EncodedChunkOutput;

    auto read_encoded_chunks(
        EncodedChunkOutput& output,
        gelex::GeneticMode mode,
        gelex::GenotypeMethod method,
        std::size_t chunk_size) -> void;

    gelex::Bed bed_;
    gelex::GenoObserver observer_;
    int64_t sample_size_{};
    int64_t num_variants_{};
};

}  // namespace gelex::genotype

#endif  // GELEX_DATA_GENOTYPE_READER_H_
