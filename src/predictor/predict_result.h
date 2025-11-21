#pragma once

#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <vector>

namespace gelex::predictor
{

/**
 * @brief Prediction result container
 *
 * Stores prediction results throughout the pipeline:
 * - Genetic values
 * - Covariate effects
 * - Total predictions
 * - Individual identifiers
 * - Metadata
 */
struct PredictResult
{
    // Individual identifiers
    std::vector<std::string> individual_ids;

    // Prediction components
    Eigen::VectorXd genetic_values;     // Genetic component predictions
    Eigen::VectorXd covariate_effects;  // Covariate component predictions
    Eigen::VectorXd
        total_predictions;  // Total predictions (genetic + covariate)

    // Metadata
    size_t num_matched_snps = 0;   // Number of SNPs successfully matched
    size_t num_reversed_snps = 0;  // Number of SNPs with reversed alleles
    size_t num_skipped_snps = 0;   // Number of SNPs skipped due to mismatch

    /**
     * @brief Check if result is valid
     */
    bool is_valid() const
    {
        return !individual_ids.empty() && genetic_values.size() > 0
               && covariate_effects.size() > 0 && total_predictions.size() > 0
               && static_cast<size_t>(genetic_values.size())
                      == individual_ids.size()
               && static_cast<size_t>(covariate_effects.size())
                      == individual_ids.size()
               && static_cast<size_t>(total_predictions.size())
                      == individual_ids.size();
    }

    /**
     * @brief Convert to individual ID -> total prediction map
     */
    std::unordered_map<std::string, double> to_map() const
    {
        std::unordered_map<std::string, double> result;
        for (size_t i = 0; i < individual_ids.size(); ++i)
        {
            result[individual_ids[i]]
                = total_predictions(static_cast<Eigen::Index>(i));
        }
        return result;
    }

    /**
     * @brief Convert to individual ID -> genetic value map
     */
    std::unordered_map<std::string, double> genetic_values_map() const
    {
        std::unordered_map<std::string, double> result;
        for (size_t i = 0; i < individual_ids.size(); ++i)
        {
            result[individual_ids[i]]
                = genetic_values(static_cast<Eigen::Index>(i));
        }
        return result;
    }

    /**
     * @brief Convert to individual ID -> covariate effect map
     */
    std::unordered_map<std::string, double> covariate_effects_map() const
    {
        std::unordered_map<std::string, double> result;
        for (size_t i = 0; i < individual_ids.size(); ++i)
        {
            result[individual_ids[i]]
                = covariate_effects(static_cast<Eigen::Index>(i));
        }
        return result;
    }

    /**
     * @brief Clear all results
     */
    void clear()
    {
        individual_ids.clear();
        genetic_values.resize(0);
        covariate_effects.resize(0);
        total_predictions.resize(0);
        num_matched_snps = 0;
        num_reversed_snps = 0;
        num_skipped_snps = 0;
    }

    /**
     * @brief Get SNP matching statistics
     */
    std::string get_matching_stats() const
    {
        return std::format(
            "Matched: {}, Reversed: {}, Skipped: {}",
            num_matched_snps,
            num_reversed_snps,
            num_skipped_snps);
    }
};

}  // namespace gelex::predictor
