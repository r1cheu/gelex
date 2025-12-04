#include "gelex/data/data_pipe.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/sample_manager.h"
#include "gelex/exception.h"
#include "loader/ccovariate_loader.h"
#include "loader/phenotype_loader.h"
#include "loader/qcovariate_loader.h"

namespace gelex
{
DataPipe::~DataPipe() = default;
DataPipe::DataPipe(const Config& config)
{
    auto fam_path = config.bed_path;
    fam_path.replace_extension(".fam");
    sample_manager_ = std::make_shared<SampleManager>(fam_path);

    if (config.phenotype_path.empty())
    {
        throw ArgumentValidationException("Phenotype file path is required.");
    }

    load_phenotype(config);

    if (!config.qcovar_path.empty())
    {
        load_qcovariates(config);
    }

    if (!config.covar_path.empty())
    {
        load_covariates(config);
    }

    intersect();
    load_additive(config);
    load_dominance(config);
    convert_to_matrices();
}

void DataPipe::load_phenotype(const Config& config)
{
    phenotype_loader_ = std::make_unique<detail::PhenotypeLoader>(
        config.phenotype_path, config.phenotype_column, config.iid_only);
    phenotype_name_ = phenotype_loader_->name();
}

void DataPipe::load_qcovariates(const Config& config)
{
    qcovar_loader_ = std::make_unique<detail::QcovarLoader>(
        config.qcovar_path, config.iid_only);
    qcovariate_names_ = qcovar_loader_->names();
}

void DataPipe::load_covariates(const Config& config)
{
    covar_loader_ = std::make_unique<detail::CCovarLoader>(
        config.covar_path, config.iid_only);
    covariate_names_ = covar_loader_->names();
}

void DataPipe::load_additive(const Config& config)
{
    load_genotype_impl<gelex::HardWenbergProcessor>(
        config, ".add", additive_matrix_);
}

void DataPipe::load_dominance(const Config& config)
{
    if (config.use_dominance_effect)
    {
        load_genotype_impl<gelex::DominantOrthogonalHWEProcessor>(
            config, ".dom", dominance_matrix_);
    }
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

    phenotype_ = phenotype_loader_ ? phenotype_loader_->load(id_map)
                                   : Eigen::VectorXd::Zero();

    Eigen::MatrixXd qcovariates;
    if (qcovar_loader_)
    {
        qcovariates = qcovar_loader_->load(id_map);
    }

    Eigen::MatrixXd covariates;
    if (covar_loader_)
    {
        covariates = covar_loader_->load(id_map);
    }
    const Eigen::Index num_fixed_effects
        = 1 + qcovariates.cols() + covariates.cols();

    fixed_effects_ = Eigen::MatrixXd::Ones(num_samples, num_fixed_effects);

    Eigen::Index current_col = 1;
    if (qcovariates.size() != 0)
    {
        fixed_effects_.middleCols(current_col, qcovariates.cols())
            = qcovariates;
        current_col += qcovariates.cols();
    }

    if (covariates.size() != 0)
    {
        fixed_effects_.middleCols(current_col, covariates.cols()) = covariates;
    }

    fixed_effect_names_.clear();
    fixed_effect_names_.emplace_back("intercept");
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
