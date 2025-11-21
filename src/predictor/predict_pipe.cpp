#include "predict_pipe.h"

#include <expected>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader.h"
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"

namespace gelex
{

auto PredictDataPipe::create(
    const Config& config,
    std::shared_ptr<SampleManager> sample_manager)
    -> std::expected<PredictDataPipe, Error>
{
    PredictDataPipe pipe;
    pipe.sample_manager_ = std::move(sample_manager);

    if (!config.qcovar_path.empty())
    {
        if (auto result = pipe.load_qcovariates(config); !result)
        {
            return std::unexpected(result.error());
        }
    }

    if (!config.covar_path.empty())
    {
        if (auto result = pipe.load_covariates(config); !result)
        {
            return std::unexpected(result.error());
        }
    }

    pipe.intersect();
    pipe.format_covariates();

    return pipe;
}

auto PredictDataPipe::load_qcovariates(const Config& config)
    -> std::expected<void, Error>
{
    auto loader
        = detail::QcovarLoader::create(config.qcovar_path, config.iid_only);

    if (!loader)
    {
        return std::unexpected(loader.error());
    }

    qcovar_loader_ = std::make_unique<detail::QcovarLoader>(std::move(*loader));
    qcovariate_names_ = qcovar_loader_->names();

    return {};
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

    if (covar_loader_)
    {
        auto key_views = create_key_views(covar_loader_->data());
        sample_manager_->intersect(key_views);
    }

    sample_manager_->finalize();
}

void PredictDataPipe::format_covariates()
{
    const auto num_samples
        = static_cast<Eigen::Index>(sample_manager_->num_common_samples());
    const auto& id_map = sample_manager_->common_id_map();

    Eigen::MatrixXd qcovariates;
    if (qcovar_loader_)
    {
        qcovariates = qcovar_loader_->load(id_map);
    }

    qcovariates_ = Eigen::MatrixXd::Ones(num_samples, qcovariates.cols() + 1);
    if (qcovariates.size() != 0)
    {
        qcovariates_.middleCols(1, qcovariates.cols()) = qcovariates;
    }
}

const std::vector<std::string>& PredictDataPipe::qcovariate_names() const
{
    return qcovariate_names_;
}

const std::vector<std::string>& PredictDataPipe::covariate_names() const
{
    return covariate_names_;
}

}  // namespace gelex
