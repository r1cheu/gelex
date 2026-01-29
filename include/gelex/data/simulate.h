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

#ifndef GELEX_DATA_SIMULATE_H_
#define GELEX_DATA_SIMULATE_H_

#include <filesystem>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>

#include <Eigen/Core>

#include "gelex/data/bed_pipe.h"
#include "gelex/data/sample_manager.h"
#include "gelex/exception.h"

namespace gelex
{

/**
 * @class PhenotypeSimulator
 * @brief Simulates phenotypes based on genetic data and specified parameters.
 *
 * This class orchestrates the simulation of quantitative traits using genotype
 * data from a PLINK BED file, a list of causal variants, and a given
 * heritability.
 */
class PhenotypeSimulator
{
   public:
    struct Config
    {
        std::filesystem::path bed_path;
        std::filesystem::path causal_variants_path;
        double heritability{0.5};
        int seed{-1};
        std::filesystem::path output_path;
    };

    explicit PhenotypeSimulator(Config config);

    PhenotypeSimulator(const PhenotypeSimulator&) = delete;
    PhenotypeSimulator(PhenotypeSimulator&&) noexcept = default;
    PhenotypeSimulator& operator=(const PhenotypeSimulator&) = delete;
    PhenotypeSimulator& operator=(PhenotypeSimulator&&) noexcept = default;
    ~PhenotypeSimulator() = default;

    /**
     * @brief Runs the phenotype simulation process.
     */
    void simulate();

   private:
    // --- Helper Methods for Simulation Stages ---
    void initialize_rng();
    auto load_or_generate_causal_effects()
        -> std::unordered_map<std::string, double>;

    static auto calculate_genetic_values(
        BedPipe& bed_pipe,
        const std::unordered_map<std::string, double>& causal_effects)
        -> Eigen::VectorXd;

    auto generate_phenotypes(const Eigen::VectorXd& genetic_values)
        -> Eigen::VectorXd;

    void write_results(
        const Eigen::VectorXd& phenotypes,
        const std::shared_ptr<SampleManager>& sample_manager) const;

    Config config_;
    std::mt19937_64 rng_;

    // A reasonable default chunk size for processing variants
    static constexpr Eigen::Index SNP_CHUNK_SIZE = 10000;
};

}  // namespace gelex

#endif  // GELEX_DATA_SIMULATE_H_
