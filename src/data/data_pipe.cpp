#include "gelex/data/data_pipe.h"

#include <algorithm>
#include <expected>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader.h"
#include "gelex/error.h"
#include "gelex/logger.h"

namespace gelex
{

auto DataPipe::create(const Config& config) -> std::expected<DataPipe, Error>
{
    DataPipe pipe;

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

    if (!config.fam_path.empty())
    {
        if (auto result = pipe.load_fam(config); !result)
        {
            return std::unexpected(result.error());
        }
    }

    pipe.intersect_ids();
    pipe.build_id_map();

    pipe.convert_to_matrices();

    return pipe;
}

auto DataPipe::load_phenotype(const Config& config)
    -> std::expected<void, Error>
{
    auto loader = detail::PhenotypeLoader::create(
        config.phenotype_path.string(),
        config.phenotype_column,
        config.iid_only);

    if (!loader)
    {
        return std::unexpected(loader.error());
    }

    phenotype_loader_
        = std::make_unique<detail::PhenotypeLoader>(std::move(*loader));
    phenotype_name_ = phenotype_loader_->phenotype_name();

    return {};
}

auto DataPipe::load_qcovariates(const Config& config)
    -> std::expected<void, Error>
{
    auto loader = detail::QcovarLoader::create(
        config.qcovar_path.string(), config.iid_only);

    if (!loader)
    {
        return std::unexpected(loader.error());
    }

    qcovar_loader_ = std::make_unique<detail::QcovarLoader>(std::move(*loader));
    qcovariate_names_ = qcovar_loader_->covariate_names();

    return {};
}

auto DataPipe::load_covariates(const Config& config)
    -> std::expected<void, Error>
{
    auto loader = detail::CovarLoader::create(
        config.covar_path.string(), config.iid_only);

    if (!loader)
    {
        return std::unexpected(loader.error());
    }

    covar_loader_ = std::make_unique<detail::CovarLoader>(std::move(*loader));
    covariate_names_ = covar_loader_->covariate_names();

    return {};
}

auto DataPipe::load_fam(const Config& config) -> std::expected<void, Error>
{
    auto loader
        = detail::FamLoader::create(config.fam_path.string(), config.iid_only);

    if (!loader)
    {
        return std::unexpected(loader.error());
    }

    fam_loader_ = std::make_unique<detail::FamLoader>(std::move(*loader));

    return {};
}

void DataPipe::intersect_ids()
{
    std::unordered_set<std::string> intersection;

    if (phenotype_loader_)
    {
        for (const auto& [id, _] : phenotype_loader_->phenotype_data())
        {
            intersection.insert(id);
        }
    }

    if (qcovar_loader_)
    {
        for (const auto& [id, _] : qcovar_loader_->covariate_data())
        {
            intersection.insert(id);
        }
    }

    if (covar_loader_)
    {
        for (const auto& [id, _] : covar_loader_->covariate_data())
        {
            intersection.insert(id);
        }
    }

    if (fam_loader_)
    {
        for (const auto& id : fam_loader_->sample_ids())
        {
            intersection.insert(id);
        }
    }

    if (phenotype_loader_)
    {
        phenotype_loader_->intersect(intersection);
    }

    if (qcovar_loader_)
    {
        qcovar_loader_->intersect(intersection);
    }

    if (covar_loader_)
    {
        covar_loader_->intersect(intersection);
    }

    if (fam_loader_)
    {
        fam_loader_->intersect(intersection);
    }

    sample_ids_.assign(intersection.begin(), intersection.end());
    std::ranges::sort(sample_ids_);

    auto logger = gelex::logging::get();
    logger->info(
        "{} common samples available for analysis after intersection",
        sample_ids_.size());
}

void DataPipe::build_id_map()
{
    id_map_.clear();
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(sample_ids_.size());
         ++i)
    {
        id_map_.emplace(sample_ids_[i], i);
    }
}

void DataPipe::convert_to_matrices()
{
    const auto num_samples = static_cast<Eigen::Index>(sample_ids_.size());

    phenotype_ = phenotype_loader_ ? std::move(phenotype_loader_->load(id_map_))
                                   : Eigen::VectorXd::Zero();

    Eigen::MatrixXd qcovariates;
    if (qcovar_loader_)
    {
        qcovariates = std::move(qcovar_loader_->load(id_map_));
    }

    Eigen::MatrixXd covariates;
    if (covar_loader_)
    {
        covariates = std::move(covar_loader_->load(id_map_));
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
