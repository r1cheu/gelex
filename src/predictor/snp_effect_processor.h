#pragma once

#include <cmath>
#include <expected>
#include <string>
#include <vector>

#include "Eigen/Core"
#include "gelex/error.h"
#include "snp_matcher.h"

namespace gelex
{

/// Structure to store column indices for .snp.eff file parsing
struct ColumnIndices
{
    int id = -1;
    int a1 = -1;
    int a2 = -1;
    int a1frq = -1;
    int add = -1;
    int dom = -1;

    /// Check if all required columns are present
    bool has_required_columns() const
    {
        return id != -1 && a1 != -1 && a2 != -1 && a1frq != -1 && add != -1;
    }
};

class SnpEffectProcessor
{
   public:
    // snp_eff_path: Path to the .snp.eff file
    static auto create(const std::filesystem::path& snp_eff_path)
        -> std::expected<SnpEffectProcessor, Error>;

    static double calculate_gevi(int genotype, const SnpEffect& info)
    {
        double p = info.p_freq;
        double q = 1.0 - p;

        // Handle monomorphic cases (p=0 or p=1)
        if (p <= 0.0 || p >= 1.0)
        {
            // For monomorphic SNPs, all genotypes are the same, so effect is
            // zero
            return 0.0;
        }

        // Additive encoding: genotype count (0, 1, 2)
        auto add_encoded = static_cast<double>(genotype);

        // Dominant encoding: based on genotype
        double dom_encoded = 0.0;
        if (genotype == 0)
        {
            dom_encoded = 0.0;
        }
        else if (genotype == 1)
        {
            dom_encoded = 2.0 * p;
        }
        else
        {  // genotype == 2
            dom_encoded = (4.0 * p) - 2.0;
        }

        // Standardized additive effect
        double add_std = (add_encoded - 2.0 * p) / std::sqrt(2.0 * p * q);

        // Standardized dominant effect
        double dom_std = (dom_encoded - (2.0 * p * p)) / (2.0 * p * q);

        // Total genetic value (GEVI)
        return (add_std * info.add_effect) + (dom_std * info.dom_effect);
    }

    /// Calculate genetic values for multiple genotypes
    /// @param genotypes Vector of genotype values (0, 1, 2)
    /// @param info SNP information
    /// @return Vector of calculated genetic values
    static std::vector<double> calculate_gevi_batch(
        const std::vector<int>& genotypes,
        const SnpEffect& info)
    {
        std::vector<double> results;
        results.reserve(genotypes.size());

        for (int genotype : genotypes)
        {
            results.push_back(calculate_gevi(genotype, info));
        }

        return results;
    }
    /// Calculate total genetic value across multiple SNPs
    /// @param genotypes Vector of genotype vectors (one per SNP)
    /// @param snp_infos Vector of SNP information
    /// @return Vector of total genetic values (one per individual)
    static std::vector<double> calculate_total_genetic_value(
        const std::vector<std::vector<int>>& genotypes,
        const std::vector<SnpEffect>& snp_infos);

    const std::vector<SnpEffect>& snp_effects() const { return snp_effects_; }

   private:
    SnpEffectProcessor() = default;
    static std::expected<std::vector<SnpEffect>, Error> parse_snp_eff_file(
        const std::string& snp_eff_file);

    static ColumnIndices assign_column_indices(
        std::span<const std::string_view> header_columns);

    static auto create_snp_info(
        std::span<const std::string_view> columns,
        const ColumnIndices& indices) -> std::expected<SnpEffect, Error>;

    std::vector<SnpEffect> snp_effects_;
    Eigen::VectorXd add_;
    Eigen::VectorXd dom_;
    bool has_dominant_effects_ = false;
};

}  // namespace gelex
