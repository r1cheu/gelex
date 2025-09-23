#include "gelex/data/data_pipe.h"

#include <expected>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader.h"
#include "gelex/error.h"

namespace gelex
{

auto DataPipe::create(
    const Config& config,
    std::shared_ptr<SampleManager> sample_manager)
    -> std::expected<DataPipe, Error>
{
    DataPipe pipe(std::move(sample_manager));

    if (!config.phenotype_path.empty())
    {
        if (auto result = pipe.load_phenotype(config); !result)
        {
            return std::unexpected(result.error());
        }
    }

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
    pipe.convert_to_matrices();

    return pipe;
}

auto DataPipe::load_phenotype(const Config& config)
    -> std::expected<void, Error>
{
    auto loader = detail::PhenotypeLoader::create(
        config.phenotype_path, config.phenotype_column, config.iid_only);

    if (!loader)
    {
        return std::unexpected(loader.error());
    }

    phenotype_loader_
        = std::make_unique<detail::PhenotypeLoader>(std::move(*loader));
    phenotype_name_ = phenotype_loader_->name();

    return {};
}

auto DataPipe::load_qcovariates(const Config& config)
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

auto DataPipe::load_covariates(const Config& config)
    -> std::expected<void, Error>
{
    auto loader
        = detail::CovarLoader::create(config.covar_path, config.iid_only);

    if (!loader)
    {
        return std::unexpected(loader.error());
    }

    covar_loader_ = std::make_unique<detail::CovarLoader>(std::move(*loader));
    covariate_names_ = covar_loader_->names();

    return {};
}

void DataPipe::intersect()
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

    if (phenotype_loader_)
    {
        auto key_views = create_key_views(phenotype_loader_->data());
        sample_manager_->intersect(key_views);
    }

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

void DataPipe::convert_to_matrices()
{
    const auto num_samples
        = static_cast<Eigen::Index>(sample_manager_->num_common_samples());
    const auto& id_map = sample_manager_->common_id_map();

    phenotype_ = phenotype_loader_ ? std::move(phenotype_loader_->load(id_map))
                                   : Eigen::VectorXd::Zero();

    Eigen::MatrixXd qcovariates;
    if (qcovar_loader_)
    {
        qcovariates = std::move(qcovar_loader_->load(id_map));
    }

    Eigen::MatrixXd covariates;
    if (covar_loader_)
    {
        covariates = std::move(covar_loader_->load(id_map));
    }
    const Eigen::Index num_fixed_effects
        = 1 + qcovariates.cols() + covariates.cols();

    fixed_effects_ = Eigen::MatrixXd::Ones(num_samples, num_fixed_effects);

    if (qcovariates.size() != 0)
    {
        fixed_effects_.middleCols(1, qcovariates.cols()) = qcovariates;
    }

    if (covariates.size() != 0)
    {
        fixed_effects_.middleCols(1 + qcovariates.cols(), covariates.cols())
            = covariates;
    }

    fixed_effect_names_.clear();
    fixed_effect_names_.emplace_back("Intercept");
    fixed_effect_names_.insert(
        fixed_effect_names_.end(),
        qcovariate_names_.begin(),
        qcovariate_names_.end());
    fixed_effect_names_.insert(
        fixed_effect_names_.end(),
        covariate_names_.begin(),
        covariate_names_.end());
}

const std::string& DataPipe::phenotype_name() const
{
    return phenotype_name_;
}

const std::vector<std::string>& DataPipe::qcovariate_names() const
{
    return qcovariate_names_;
}

const std::vector<std::string>& DataPipe::covariate_names() const
{
    return covariate_names_;
}

const std::vector<std::string>& DataPipe::fixed_effect_names() const
{
    return fixed_effect_names_;
}

}  // namespace gelex
