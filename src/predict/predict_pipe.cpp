#include "predict_pipe.h"

#include <filesystem>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader/qcovariate_loader.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/sample_manager.h"

namespace gelex
{

PredictDataPipe::PredictDataPipe(const Config& config)
{
    auto fam_path = config.bed_path;
    fam_path.replace_extension(".fam");
    sample_manager_
        = std::make_shared<SampleManager>(fam_path, config.iid_only);

    if (!config.qcovar_path.empty())
    {
        load_qcovariates(config);
    }

    if (!config.dcovar_path.empty())
    {
        load_dcovariates(config);
    }

    intersect();
    format_dcovariates();

    load_genotype(config);
}

void PredictDataPipe::load_qcovariates(const Config& config)
{
    qcovar_loader_ = std::make_unique<detail::QuantitativeCovariateLoader>(
        config.qcovar_path, config.iid_only);
    qcovariate_names_ = qcovar_loader_->names();
}

void PredictDataPipe::load_dcovariates(const Config& config)
{
    dcovar_loader_ = std::make_unique<detail::DcovarPredictLoader>(
        config.dcovar_path, config.iid_only);
    dcovariate_names_ = dcovar_loader_->names();
}

void PredictDataPipe::load_genotype(const Config& config)
{
    BedPipe bed_pipe(config.bed_path, sample_manager_);
    genotypes_ = bed_pipe.load();
}

void PredictDataPipe::intersect()
{
    auto create_key_views = [](const auto& data)
    {
        std::vector<std::string_view> key_views;
        key_views.reserve(data.size());
        for (const auto& [k, v] : data)
        {
            key_views.emplace_back(k);
        }
        return key_views;
    };

    if (qcovar_loader_)
    {
        auto key_views = create_key_views(qcovar_loader_->data());
        sample_manager_->intersect(key_views);
    }

    if (dcovar_loader_)
    {
        auto key_views = create_key_views(dcovar_loader_->data());
        sample_manager_->intersect(key_views);
    }

    sample_manager_->finalize();
}

void PredictDataPipe::format_dcovariates()
{
    const auto num_samples
        = static_cast<Eigen::Index>(sample_manager_->num_common_samples());
    const auto& id_map = sample_manager_->common_id_map();
    sample_ids_ = sample_manager_->common_ids();

    if (qcovar_loader_)
    {
        auto qcov = qcovar_loader_->load(id_map);
        qcovariates_ = Eigen::MatrixXd::Ones(num_samples, qcov.X.cols() + 1);
        if (qcov.X.size() != 0)
        {
            qcovariates_.middleCols(1, qcov.X.cols()) = qcov.X;
        }
    }
    else
    {
        qcovariates_ = Eigen::MatrixXd::Ones(num_samples, 1);
    }

    if (dcovar_loader_)
    {
        dcovariates_ = dcovar_loader_->load(id_map);
    }
}

}  // namespace gelex
