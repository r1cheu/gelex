#pragma once

#include <expected>
#include <filesystem>
#include <vector>

#include "Eigen/Core"
#include "covariate_processor.h"
#include "gelex/error.h"
#include "predict_result.h"
#include "snp_effect_processor.h"
#include "snp_matcher.h"

namespace gelex
{

class PredictorEngine
{
   public:
    static auto create(
        const std::filesystem::path& snp_eff_path,
        const std::filesystem::path& covariate_param_path)
        -> std::expected<PredictorEngine, Error>;

    auto predict(
        const std::filesystem::path& prediction_bed_prefix,
        const std::vector<detail::IndividualData>& covariate_data = {})
        -> std::expected<predictor::PredictResult, Error>;

    const std::vector<SnpEffect>& snp_effects() const
    {
        return snp_effect_processor_.snp_effects();
    }
    const MatchInfo& match_info() const { return match_info_; }
    const Eigen::MatrixXd& genotypes() const { return genotypes_; }

   private:
    PredictorEngine(
        SnpEffectProcessor snp_effect_processor,
        detail::CovariateProcessor covariate_processor)
        : snp_effect_processor_(std::move(snp_effect_processor)),
          covariate_processor_(std::move(covariate_processor))
    {
    }

    SnpEffectProcessor snp_effect_processor_;
    detail::CovariateProcessor covariate_processor_;

    MatchInfo match_info_;
    Eigen::MatrixXd genotypes_;
    predictor::PredictResult last_result_;
};

}  // namespace gelex
