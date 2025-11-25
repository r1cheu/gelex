#pragma once

#include <filesystem>
#include <fstream>

#include <Eigen/Core>

#include "../src/data/loader.h"
#include "gelex/types/mcmc_results.h"

namespace gelex
{
namespace detail
{

class BimLoader;
}

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
