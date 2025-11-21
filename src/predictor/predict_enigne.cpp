#include "predict_enigne.h"

#include <expected>
#include <filesystem>

#include "snp_effect_processor.h"

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
}

}  // namespace gelex
