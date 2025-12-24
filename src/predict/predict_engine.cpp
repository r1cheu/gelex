#include "predict_engine.h"

#include <filesystem>
#include <memory>

#include "gelex/exception.h"
#include "gelex/logger.h"
#include "predict/genotype_aligner.h"
#include "predict_params_pipe.h"
#include "predict_pipe.h"

namespace gelex
{

void PredictEngine::Config::validate() const
{
    if (bed_path.empty())
    {
        throw InvalidInputException("BED path must be provided");
    }
    if (snp_effect_path.empty())
    {
        throw InvalidInputException("SNP effect path must be provided");
    }
    if (covar_effect_path.empty())
    {
        throw InvalidInputException("Covariate effect path must be provided");
    }
}

PredictEngine::PredictEngine(const Config& config) : config_(config)
{
    config_.validate();

    load_parameters();
    load_data();
    validate_dimensions();
}

void PredictEngine::load_parameters()
{
    PredictParamsPipe::Config params_config{
        .snp_effect_path = config_.snp_effect_path,
        .covar_effect_path = config_.covar_effect_path};

    PredictParamsPipe params_pipe(params_config);

    snp_effects_ = std::move(params_pipe).take_snp_effects();
    covar_effects_ = std::move(params_pipe).take_covar_effects();
}

void PredictEngine::load_data()
{
    PredictDataPipe::Config data_config{
        .bed_path = config_.bed_path,
        .qcovar_path = config_.qcovar_path,
        .dcovar_path = config_.dcovar_path,
        .iid_only = config_.iid_only};

    PredictDataPipe data_pipe(data_config);
    data_ = std::move(data_pipe).take_data();
    GenotypeAligner genotype_filter(config_.bed_path, snp_effects_);
    data_.genotype = genotype_filter.load(std::move(data_.genotype));
}

void PredictEngine::run() {}

void PredictEngine::compute_predictions() {}

void PredictEngine::write_output() {}

}  // namespace gelex
