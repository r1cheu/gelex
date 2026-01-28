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

#ifndef GELEX_ESTIMATOR_BAYES_SNP_QUANT_GENETIC_WRITER_H_
#define GELEX_ESTIMATOR_BAYES_SNP_QUANT_GENETIC_WRITER_H_

#include <filesystem>
#include <fstream>

#include <Eigen/Core>

#include "../src/data/loader/bim_loader.h"
#include "gelex/types/mcmc_results.h"

namespace gelex
{

class SnpQuantGeneticWriter
{
   public:
    /**
     * @brief Constructor for SnpQuantGeneticWriter
     *
     * @param result MCMC result containing SNP effects
     * @param bim_file_path Path to BIM file for SNP information
     */
    explicit SnpQuantGeneticWriter(
        const MCMCResult& result,
        const std::filesystem::path& bim_file_path);

    /**
     * @brief Write original genetic effects to file
     *
     * @param path Output file path
     */
    void write(const std::filesystem::path& path) const;

   private:
    const MCMCResult* result_;
    detail::BimLoader bim_loader_;

    // Helper methods for writing different SNP effect components
    void write_header(std::ofstream& stream) const;
    void write_snp_row(std::ofstream& stream, Eigen::Index snp_index) const;
    double write_snp_basic_info(std::ofstream& stream, Eigen::Index snp_index)
        const;
};

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_SNP_QUANT_GENETIC_WRITER_H_
