#include "predict_engine.h"

#include <format>

#include "covar_predictor.h"
#include "gelex/exception.h"
#include "predict/genotype_aligner.h"
#include "predict_params_pipe.h"
#include "predict_pipe.h"
#include "predict_writer.h"
#include "snp_predictor.h"

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
    if (output_path.empty())
    {
        throw InvalidInputException("Output path must be provided");
    }
}

PredictEngine::PredictEngine(const Config& config) : config_(config)
{
    config_.validate();
}

void PredictEngine::run()
{
    load_parameters();
    load_data();
    validate_dimensions();
    compute();
    write();
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
    sample_ids_ = data_.sample_ids;

    GenotypeAligner genotype_filter(config_.bed_path, snp_effects_);
    data_.genotype = genotype_filter.align(std::move(data_.genotype));
}

void PredictEngine::validate_dimensions()
{
    const Eigen::Index n_snps = data_.genotype.cols();
    const auto n_snp_effects = static_cast<Eigen::Index>(snp_effects_.size());

    if (n_snps != n_snp_effects)
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: genotype matrix has {} SNPs, but SNP "
                "effects has {}",
                n_snps,
                n_snp_effects));
    }
}

void PredictEngine::compute()
{
    SnpPredictor snp_predictor(snp_effects_);
    auto snp_result = snp_predictor.compute(data_.genotype);
    snp_predictions_ = snp_result.total();
    add_predictions_ = std::move(snp_result.add);
    dom_predictions_ = std::move(snp_result.dom);

    CovarPredictor covar_predictor(covar_effects_);
    auto covar_result = covar_predictor.compute(data_);
    covar_predictions_ = std::move(covar_result.predictions);
    covar_prediction_names_ = std::move(covar_result.names);
    predictions_ = snp_predictions_ + covar_predictions_.rowwise().sum();
}

void PredictEngine::write()
{
    PredictWriter writer(config_.output_path, config_.iid_only);
    writer.write(
        predictions_,
        sample_ids_,
        add_predictions_,
        dom_predictions_,
        covar_predictions_,
        covar_prediction_names_);
}

}  // namespace gelex
