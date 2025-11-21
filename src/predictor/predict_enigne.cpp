#include "predictor_enigne.h"

#include <expected>
#include <filesystem>

#include "Eigen/Core"
#include "genotype_reader.h"
#include "snp_effect_processor.h"
#include "snp_matcher.h"

namespace gelex
{

auto PredictorEngine::create(
    const std::filesystem::path& snp_eff_path,
    const std::filesystem::path& covariate_param_path)
    -> std::expected<PredictorEngine, Error>
{
    // Load SNP effects
    auto snp_effect_result = SnpEffectProcessor::create(snp_eff_path);
    if (!snp_effect_result)
    {
        return std::unexpected(snp_effect_result.error());
    }

    // Load covariate parameters
    auto covariate_result
        = detail::CovariateProcessor::create(covariate_param_path);
    if (!covariate_result)
    {
        return std::unexpected(covariate_result.error());
    }

    return PredictorEngine(
        std::move(snp_effect_result.value()),
        std::move(covariate_result.value()));
}

auto PredictorEngine::predict(
    const std::filesystem::path& prediction_bed_prefix,
    const std::vector<detail::IndividualData>& covariate_data)
    -> std::expected<predictor::PredictResult, Error>
{
    // Create genotype reader
    auto genotype_reader_result
        = PredictorGenotypeReader::create(prediction_bed_prefix);
    if (!genotype_reader_result)
    {
        return std::unexpected(genotype_reader_result.error());
    }
    auto& genotype_reader = genotype_reader_result.value();

    // Match SNPs between effect file and prediction genotypes
    auto snp_matcher_result = SnpMatcher::create(prediction_bed_prefix);
    if (!snp_matcher_result)
    {
        return std::unexpected(snp_matcher_result.error());
    }
    auto& snp_matcher = snp_matcher_result.value();

    match_info_ = snp_matcher.match(snp_effect_processor_.snp_effects());

    // Calculate SNP matching statistics
    size_t num_matched = 0;
    size_t num_reversed = 0;
    size_t num_skipped = 0;
    for (const auto& match_type : match_info_.read_plan)
    {
        switch (match_type)
        {
            case MatchType::keep:
                ++num_matched;
                break;
            case MatchType::reverse:
                ++num_reversed;
                break;
            case MatchType::skip:
                ++num_skipped;
                break;
        }
    }

    // Read and process genotypes
    auto genotype_result = genotype_reader.process(match_info_);
    if (!genotype_result)
    {
        return std::unexpected(genotype_result.error());
    }
    genotypes_ = std::move(genotype_result.value());

    // Calculate genetic values: genotypes * matched SNP effects
    Eigen::VectorXd genetic_values = genotypes_;

    // Calculate covariate effects
    Eigen::VectorXd covariate_effects(genotype_reader.num_samples());
    if (covariate_data.empty())
    {
        // Use default covariate prediction (just intercept)
        for (Eigen::Index i = 0; i < covariate_effects.size(); ++i)
        {
            covariate_effects[i]
                = covariate_processor_.predict(detail::IndividualData{});
        }
    }
    else
    {
        // Use provided covariate data
        if (covariate_data.size()
            != static_cast<size_t>(genotype_reader.num_samples()))
        {
            return std::unexpected(
                Error{
                    ErrorCode::InvalidData,
                    "Number of covariate data entries does not match number of "
                    "samples"});
        }
        for (size_t i = 0; i < covariate_data.size(); ++i)
        {
            covariate_effects[static_cast<Eigen::Index>(i)]
                = covariate_processor_.predict(covariate_data[i]);
        }
    }

    // Calculate total predictions
    Eigen::VectorXd total_predictions = genetic_values + covariate_effects;

    // Get individual IDs (placeholder - would need to read from FAM file)
    std::vector<std::string> individual_ids(genotype_reader.num_samples());
    for (size_t i = 0; i < individual_ids.size(); ++i)
    {
        individual_ids[i] = "sample_" + std::to_string(i);
    }

    // Store and return results
    last_result_ = predictor::PredictResult{
        std::move(individual_ids),
        std::move(genetic_values),
        std::move(covariate_effects),
        std::move(total_predictions),
        num_matched,
        num_reversed,
        num_skipped};

    return last_result_;
}

}  // namespace gelex
