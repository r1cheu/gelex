#include "gelex/data/data_pipe.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "../src/estimator/bayes/indicator.h"
#include "barkeep.h"
#include "gelex/data/sample_manager.h"
#include "gelex/data/variant_processor.h"
#include "gelex/exception.h"
#include "loader/dcovariate_loader.h"
#include "loader/phenotype_loader.h"
#include "loader/qcovariate_loader.h"

namespace bk = barkeep;

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
    num_genotype_samples_ = sample_manager_->num_common_samples();
}

PhenoStats DataPipe::load_phenotypes()
{
    if (config_.phenotype_path.empty())
    {
        throw ArgumentValidationException("Phenotype file path is required.");
    }

    phenotype_loader_ = std::make_unique<detail::PhenotypeLoader>(
        config_.phenotype_path, config_.phenotype_column, config_.iid_only);

    return PhenoStats{
        .samples_loaded = phenotype_loader_->data().size(),
        .trait_name = phenotype_loader_->name()};
}

CovarStats DataPipe::load_covariates()
{
    std::vector<std::string> q_names;
    std::vector<std::string> d_names;

    if (!config_.qcovar_path.empty())
    {
        qcovar_loader_ = std::make_unique<detail::QcovarLoader>(
            config_.qcovar_path, config_.iid_only);
        q_names = qcovar_loader_->names();
    }

    if (!config_.dcovar_path.empty())
    {
        dcovar_loader_ = std::make_unique<detail::DcovarLoader>(
            config_.dcovar_path, config_.iid_only);
        d_names = dcovar_loader_->names();
    }

    return CovarStats{
        .qcovar_loaded = q_names.size(),
        .dcovar_loaded = d_names.size(),
        .q_names = std::move(q_names),
        .d_names = std::move(d_names)};
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
    auto intersect = [&](const auto& data)
    {
        total_before = std::max(total_before, data.size());
        auto key_views = create_key_views(data);
        sample_manager_->intersect(key_views);
    };

    if (phenotype_loader_)
    {
        intersect(phenotype_loader_->data());
    }

    if (qcovar_loader_)
    {
        intersect(qcovar_loader_->data());
    }

    if (dcovar_loader_)
    {
        intersect(dcovar_loader_->data());
    }

    sample_manager_->finalize();

    size_t common = sample_manager_->num_common_samples();

    return IntersectionStats{
        .total_samples = total_before,
        .common_samples = common,
        .excluded_samples = total_before - common};
}

GenotypeStats DataPipe::load_additive_matrix()
{
    load_genotype_impl<gelex::HardWenbergProcessor>(".add", additive_matrix_);

    int64_t num_mono_snps = 0;
    int64_t total_snps = 0;

    if (additive_matrix_)
    {
        std::visit(
            [&](auto&& arg)
            {
                num_mono_snps = arg.num_mono();
                total_snps = arg.cols();
            },
            *additive_matrix_);
    }

    return GenotypeStats{
        .num_snps = total_snps, .monomorphic_snps = num_mono_snps};
}

GenotypeStats DataPipe::load_dominance_matrix()
{
    load_genotype_impl<gelex::DominantOrthogonalHWEProcessor>(
        ".dom", dominance_matrix_);

    int64_t num_mono_snps = 0;
    int64_t total_snps = 0;

    if (dominance_matrix_)
    {
        std::visit(
            [&](auto&& arg)
            {
                num_mono_snps = arg.num_mono();
                total_snps = arg.cols();
            },
            *dominance_matrix_);
    }

    return GenotypeStats{
        .num_snps = total_snps, .monomorphic_snps = num_mono_snps};
}

void DataPipe::finalize()
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
    if (qcovar_loader_)
    {
        auto names = qcovar_loader_->names();
        fixed_effect_names_.insert(
            fixed_effect_names_.end(), names.begin(), names.end());
    }
    if (dcovar_loader_)
    {
        auto names = dcovar_loader_->names();
        fixed_effect_names_.insert(
            fixed_effect_names_.end(), names.begin(), names.end());
    }
}

const std::vector<std::string>& DataPipe::fixed_effect_names() const
{
    return fixed_effect_names_;
}

}  // namespace gelex
