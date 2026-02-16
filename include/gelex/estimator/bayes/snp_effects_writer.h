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
#include <fstream>

#include <Eigen/Core>

#include "gelex/data/genotype/loader/bim_loader.h"
#include "gelex/types/mcmc_results.h"

namespace gelex
{
class SnpEffectsWriter
{
   public:
    explicit SnpEffectsWriter(
        const MCMCResult& result,
        const std::filesystem::path& bim_file_path);

    void write(const std::filesystem::path& path) const;

   private:
    const MCMCResult* result_;
    detail::BimLoader bim_loader_;

    void write_header(std::ofstream& stream) const;
    void write_snp_row(std::ofstream& stream, Eigen::Index snp_index) const;
    void write_snp_basic_info(std::ofstream& stream, Eigen::Index snp_index)
        const;
    void write_additive_effects(std::ofstream& stream, Eigen::Index snp_index)
        const;
    void write_add_component_probabilities(
        std::ofstream& stream,
        Eigen::Index snp_index) const;
    void write_add_pip(std::ofstream& stream, Eigen::Index snp_index) const;
    void write_dom_component_probabilities(
        std::ofstream& stream,
        Eigen::Index snp_index) const;
    void write_dom_pip(std::ofstream& stream, Eigen::Index snp_index) const;
    void write_dominant_effects(std::ofstream& stream, Eigen::Index snp_index)
        const;
};

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_SNP_EFFECTS_WRITER_H_
