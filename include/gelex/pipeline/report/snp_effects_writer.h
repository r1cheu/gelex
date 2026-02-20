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

#ifndef GELEX_ESTIMATOR_BAYES_SNP_EFFECTS_WRITER_H_
#define GELEX_ESTIMATOR_BAYES_SNP_EFFECTS_WRITER_H_

#include <filesystem>
#include <memory>
#include <string>

#include <Eigen/Core>

#include "gelex/data/loader/bim_loader.h"
#include "gelex/types/mcmc_results.h"

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
        const MCMCResult& result,
        const std::filesystem::path& bim_file_path,
        const std::filesystem::path& output_path);
    ~SnpEffectsWriter();

    auto write() -> void;

   private:
    const MCMCResult* result_;
    detail::BimLoader bim_loader_;
    std::unique_ptr<detail::TextWriter> writer_;
    std::string row_buf_;

    auto write_header() -> void;
    auto write_snp_row(Eigen::Index snp_index) -> void;
    auto write_snp_basic_info(Eigen::Index snp_index) -> void;
    auto write_effects(const BaseMarkerSummary* effect, Eigen::Index snp_index)
        -> void;
    auto write_component_probs(
        const BaseMarkerSummary* effect,
        Eigen::Index snp_index) -> void;
    auto write_pip(const BaseMarkerSummary* effect, Eigen::Index snp_index)
        -> void;
};

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_SNP_EFFECTS_WRITER_H_
