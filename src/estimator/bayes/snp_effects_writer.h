#pragma once

#include <filesystem>
#include <fstream>

#include <Eigen/Core>

#include "gelex/types/mcmc_results.h"

#include "../src/data/loader/bim_loader.h"

namespace gelex
{
class SnpEffectsWriter
{
   public:
    /**
     * @brief Constructor for SnpEffectsWriter
     *
     * @param result MCMC result containing SNP effects
     * @param bim_file_path Path to BIM file for SNP information
     */
    explicit SnpEffectsWriter(
        const MCMCResult& result,
        const std::filesystem::path& bim_file_path);

    /**
     * @brief Write SNP effects to file
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
