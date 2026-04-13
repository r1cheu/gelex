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

#ifndef GELEX_MODEL_BAYES_WRITER_SNP_EFFECTS_WRITER_H_
#define GELEX_MODEL_BAYES_WRITER_SNP_EFFECTS_WRITER_H_

#include <cstdint>
#include <filesystem>
#include <memory>
#include <span>
#include <string>

#include <Eigen/Core>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/types/mcmc_result.h"

namespace gelex::detail
{
class TextWriter;
}

namespace gelex
{

class SnpEffectsWriter
{
   public:
    SnpEffectsWriter(
        const mcmc::Result& result,
        const std::filesystem::path& bim_file_path,
        const std::filesystem::path& output_path);
    ~SnpEffectsWriter();
    SnpEffectsWriter(const SnpEffectsWriter&) = delete;
    SnpEffectsWriter& operator=(const SnpEffectsWriter&) = delete;
    SnpEffectsWriter(SnpEffectsWriter&&) = delete;
    SnpEffectsWriter& operator=(SnpEffectsWriter&&) = delete;

    auto write() -> void;

   private:
    const mcmc::Result* result_;
    const GeneticSummary* additive_{};
    const GeneticSummary* dominant_{};
    df::DataFrame<std::string> bim_;
    std::unique_ptr<detail::TextWriter> writer_;
    std::string row_buf_;

    auto write_header() -> void;
    auto write_snp_row(Eigen::Index snp_index) -> void;
    auto write_snp_basic_info(Eigen::Index snp_index) -> void;
    auto write_effects(const GeneticSummary* effect, Eigen::Index snp_index)
        -> void;
    auto write_component_probs(
        const GeneticSummary* effect,
        Eigen::Index snp_index) -> void;
    auto write_pip(const GeneticSummary* effect, Eigen::Index snp_index)
        -> void;

    std::span<const std::string> bim_keys_;
    std::span<const std::string> bim_chrom_;
    std::span<const std::int32_t> bim_pos_;
    std::span<const std::string> bim_a1_;
    std::span<const std::string> bim_a2_;
};

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_WRITER_SNP_EFFECTS_WRITER_H_
