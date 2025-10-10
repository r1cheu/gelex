#pragma once

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>

#include <Eigen/Core>

#include "../src/data/snp_info_loader.h"
#include "gelex/estimator/bayes/result.h"

namespace gelex
{

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

   private:
    const MCMCResult& result_;
    std::unique_ptr<SnpInfoLoader> snp_info_loader_;

    // Helper methods for writing different file types
    void write_parameter_file(const std::filesystem::path& path) const;
    void write_snp_effects(const std::filesystem::path& path) const;

    // Helper methods for writing summary statistics
    void write_summary_statistics(
        std::ofstream& stream,
        const PosteriorSummary& stats,
        Eigen::Index n_params) const;

    // Format HPDI labels
    std::string get_hpdi_low_label() const;
    std::string get_hpdi_high_label() const;
};

}  // namespace gelex
