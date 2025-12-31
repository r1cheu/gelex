#include "gelex/data/data_pipe.h"

#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/sample_manager.h"
#include "gelex/exception.h"
#include "loader/dcovariate_loader.h"
#include "loader/phenotype_loader.h"
#include "loader/qcovariate_loader.h"

namespace gelex
{
DataPipe::~DataPipe() = default;
DataPipe::DataPipe(DataPipe&&) noexcept = default;
DataPipe& DataPipe::operator=(DataPipe&&) noexcept = default;

DataPipe::DataPipe(const Config& config) : config_(config)
{
    auto fam_path = config.bed_path;
    fam_path.replace_extension(".fam");
    sample_manager_ = std::make_shared<SampleManager>(fam_path);
}

PhenoStats DataPipe::load_phenotypes()
{
    if (config_.phenotype_path.empty())
    {
        throw ArgumentValidationException("Phenotype file path is required.");
    }

    phenotype_loader_ = std::make_unique<detail::PhenotypeLoader>(
        config_.phenotype_path, config_.phenotype_column, config_.iid_only);
    phenotype_name_ = phenotype_loader_->name();

    return PhenoStats{
        .samples_loaded = phenotype_loader_->data().size(),
        .trait_name = phenotype_name_};
}

CovarStats DataPipe::load_covariates()
{
    if (!config_.qcovar_path.empty())
    {
        qcovar_loader_ = std::make_unique<detail::QcovarLoader>(
            config_.qcovar_path, config_.iid_only);
        qcovariate_names_ = qcovar_loader_->names();
    }

    if (!config_.dcovar_path.empty())
    {
        dcovar_loader_ = std::make_unique<detail::DcovarLoader>(
            config_.dcovar_path, config_.iid_only);
        dcovariate_names_ = dcovar_loader_->names();
    }

    return CovarStats{
        .qcovar_loaded = qcovariate_names_.size(),
        .dcovar_loaded = dcovariate_names_.size(),
        .q_names = qcovariate_names_,
        .d_names = dcovariate_names_};
}

IntersectionStats DataPipe::intersect_samples()
{
    size_t total_before = sample_manager_->num_common_samples();

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

    if (dcovar_loader_)
    {
        auto key_views = create_key_views(dcovar_loader_->data());
        sample_manager_->intersect(key_views);
    }
    sample_manager_->finalize();

    size_t common = sample_manager_->num_common_samples();

    return IntersectionStats{
        .total_samples = total_before,
        .common_samples = common,
        .excluded_samples = total_before - common};
}

GenotypeStats DataPipe::build_matrices(
    std::function<void(size_t, size_t)> progress_callback)
{
    load_genotype_impl<gelex::HardWenbergProcessor>(
        ".add", additive_matrix_, progress_callback);

    if (config_.use_dominance_effect)
    {
        // For the second pass, if we want to show progress again, we might need
        // to reset or handle the callback appropriately. For now, we just pass
        // it again, so the progress bar (if stateful) might need to be reset
        // externally or the UI should handle two distinct progress phases.
        load_genotype_impl<gelex::DominantOrthogonalHWEProcessor>(
            ".dom", dominance_matrix_, progress_callback);
    }

    convert_to_matrices();

    size_t snps = 0;
    if (additive_matrix_)
    {
        std::visit([&](auto&& arg) { snps = arg.cols(); }, *additive_matrix_);
    }

    return GenotypeStats{
        .snps_loaded = snps,
        .snps_total = snps,
        .dominance_loaded = config_.use_dominance_effect};
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
    if (dcovar_loader_)
    {
        covariates = dcovar_loader_->load(id_map);
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
        dcovariate_names_.begin(),
        dcovariate_names_.end());
}

const std::string& DataPipe::phenotype_name() const
{
    return phenotype_name_;
}

const std::vector<std::string>& DataPipe::qcovariate_names() const
{
    return qcovariate_names_;
}

const std::vector<std::string>& DataPipe::dcovariate_names() const
{
    return dcovariate_names_;
}

const std::vector<std::string>& DataPipe::fixed_effect_names() const
{
    return fixed_effect_names_;
}

}  // namespace gelex
