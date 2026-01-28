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

#ifndef GELEX_ESTIMATOR_BAYES_RESULT_WRITER_H_
#define GELEX_ESTIMATOR_BAYES_RESULT_WRITER_H_

#include <filesystem>

#include "../src/estimator/bayes/parameter_writer.h"
#include "../src/estimator/bayes/snp_effects_writer.h"
#include "../src/estimator/bayes/snp_quant_genetic_writer.h"
#include "gelex/types/mcmc_results.h"

namespace gelex
{

/**
 * @brief Main result writer that coordinates specialized writers
 *
 * This class provides a facade interface that delegates to specialized
 * writers for different types of MCMC results:
 * - ParameterWriter: Handles fixed effects, random effects, variances
 * - SnpEffectsWriter: Handles SNP-specific effects with metadata
 */
class MCMCResultWriter
{
   public:
    /**
     * @brief Constructor for MCMCResultWriter
     *
     * @param result MCMC result to write
     * @param bim_file_path Path to BIM file for SNP information
     */
    explicit MCMCResultWriter(
        const MCMCResult& result,
        const std::filesystem::path& bim_file_path);

    /**
     * @brief Save MCMC results to files with the given prefix
     *
     * @param prefix Output file prefix (e.g., from CLI --out option)
     */
    void save(const std::filesystem::path& prefix) const;

    /**
     * @brief Get access to the parameter writer
     *
     * @return const ParameterWriter&
     */
    const ParameterWriter& parameter_writer() const
    {
        return parameter_writer_;
    }

    /**
     * @brief Get access to the SNP effects writer
     *
     * @return const SnpEffectsWriter&
     */
    const SnpEffectsWriter& snp_effects_writer() const
    {
        return snp_effects_writer_;
    }

    /**
     * @brief Get access to the quantitative genetic effects writer
     *
     * @return const SnpQuantGeneticWriter&
     */
    const SnpQuantGeneticWriter& snp_quant_genetic_writer() const
    {
        return snp_quant_genetic_writer_;
    }

   private:
    ParameterWriter parameter_writer_;
    SnpEffectsWriter snp_effects_writer_;
    SnpQuantGeneticWriter snp_quant_genetic_writer_;
};

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_RESULT_WRITER_H_
