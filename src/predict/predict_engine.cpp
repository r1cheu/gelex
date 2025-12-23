#include "predict_engine.h"

#include <filesystem>
#include <memory>

#include "gelex/exception.h"
#include "gelex/logger.h"
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

    auto logger = logging::get();
    if (logger)
    {
        logger->info("Initializing PredictEngine with:");
        logger->info("  BED: {}", config_.bed_path.string());
        logger->info("  SNP effects: {}", config_.snp_effect_path.string());
        logger->info(
            "  Covariate effects: {}", config_.covar_effect_path.string());
        if (!config_.qcovar_path.empty())
        {
            logger->info(
                "  Quantitative covariates: {}", config_.qcovar_path.string());
        }
        if (!config_.dcovar_path.empty())
        {
            logger->info(
                "  Discrete covariates: {}", config_.dcovar_path.string());
        }
        logger->info(
            "  IID-only mode: {}", config_.iid_only ? "true" : "false");
    }

    load_parameters();
    load_data();
    validate_dimensions();
}

void PredictEngine::load_parameters()
{
    PredictParamsPipe::Config params_config{
        .snp_effect_path = config_.snp_effect_path,
        .covar_effect_path = config_.covar_effect_path};
    params_pipe_ = std::make_unique<PredictParamsPipe>(params_config);

    snp_effects_ = std::move(*params_pipe_).take_snp_effects();
    covar_effects_ = std::move(*params_pipe_).take_covar_effects();
}

void PredictEngine::load_data()
{
    PredictDataPipe::Config data_config{
        .bed_path = config_.bed_path,
        .qcovar_path = config_.qcovar_path,
        .dcovar_path = config_.dcovar_path,
        .iid_only = config_.iid_only};

    data_pipe_ = std::make_unique<PredictDataPipe>(data_config);
}

void PredictEngine::validate_dimensions()
{
    auto logger = logging::get();
    if (logger)
    {
        const auto& snp_effects = params_pipe_->snp_effects();
        logger->info("Validation: {} SNP effects to match", snp_effects.size());
    }
}

void PredictEngine::run()
{
    auto logger = logging::get();
    if (logger)
    {
        logger->info("Starting prediction computation");
    }

    compute_predictions();
    write_output();

    if (logger)
    {
        logger->info("Prediction completed successfully");
    }
}

void PredictEngine::compute_predictions()
{
    throw InvalidOperationException("compute_predictions not yet implemented");
}

void PredictEngine::write_output()
{
    throw InvalidOperationException("write_output not yet implemented");
}

}  // namespace gelex
