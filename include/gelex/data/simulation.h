#ifndef GELEX_DATA_SIMULATION_H_
#define GELEX_DATA_SIMULATION_H_

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

#endif  // GELEX_DATA_SIMULATION_H_
