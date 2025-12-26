#include "predict_engine.h"

#include <format>
#include <fstream>
#include <ranges>

#include "../src/data/parser.h"
#include "covar_predictor.h"
#include "gelex/exception.h"
#include "predict/genotype_aligner.h"
#include "predict_params_pipe.h"
#include "predict_pipe.h"
#include "predict_writer.h"
#include "snp_predictor.h"

namespace gelex
{

namespace
{

auto load_fam_fid_iid(const std::filesystem::path& fam_path)
    -> std::pair<std::vector<std::string>, std::vector<std::string>>
{
    std::vector<std::string> fids;
    std::vector<std::string> iids;
    fids.reserve(1024);
    iids.reserve(1024);

    auto file = detail::open_file<std::ifstream>(fam_path, std::ios::in);
    std::string line;

    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        auto parts = line | std::ranges::views::split(' ')
                     | std::ranges::views::take(2);
        auto it = parts.begin();

        std::string_view fid(*it);
        fids.emplace_back(fid);

        ++it;
        std::string_view iid(*it);
        iids.emplace_back(iid);
    }

    return {std::move(fids), std::move(iids)};
}

}  // namespace

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

    auto fam_path = config_.bed_path;
    fam_path.replace_extension(".fam");
    auto [fids, iids] = load_fam_fid_iid(fam_path);

    const auto n_samples = static_cast<Eigen::Index>(sample_ids_.size());
    fids_.resize(n_samples);
    iids_.resize(n_samples);

    std::unordered_map<std::string, Eigen::Index> fid_iid_to_idx;
    fid_iid_to_idx.reserve(n_samples);

    for (Eigen::Index i = 0; i < n_samples; ++i)
    {
        fid_iid_to_idx[sample_ids_[i]] = i;
    }

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(fids.size()); ++i)
    {
        const auto& fid = fids[i];
        const auto& iid = iids[i];
        const auto key
            = config_.iid_only ? iid : std::format("{}_{}", fid, iid);

        auto it = fid_iid_to_idx.find(key);
        if (it != fid_iid_to_idx.end())
        {
            fids_[it->second] = fid;
            iids_[it->second] = iid;
        }
    }

    GenotypeAligner genotype_filter(config_.bed_path, snp_effects_);
    data_.genotype = genotype_filter.load(std::move(data_.genotype));
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
    auto snp_result = snp_predictor.compute_split(data_.genotype);
    add_predictions_ = std::move(snp_result.add);
    dom_predictions_ = std::move(snp_result.dom);
    snp_predictions_ = add_predictions_ + dom_predictions_;

    CovarPredictor covar_predictor(covar_effects_);
    auto covar_result = covar_predictor.compute(data_);
    covar_predictions_ = std::move(covar_result.predictions);
    covar_prediction_names_ = std::move(covar_result.names);

    predictions_ = snp_predictions_ + covar_predictions_.rowwise().sum();
}

void PredictEngine::write()
{
    PredictWriter writer(config_.output_path);
    writer.write(
        predictions_,
        fids_,
        iids_,
        add_predictions_,
        dom_predictions_,
        covar_predictions_,
        covar_prediction_names_);
}

}  // namespace gelex
