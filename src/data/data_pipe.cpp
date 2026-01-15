#include "gelex/data/data_pipe.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/sample_manager.h"
#include "gelex/data/variant_processor.h"
#include "gelex/exception.h"
#include "grm_loader.h"
#include "loader/dcovariate_loader.h"
#include "loader/phenotype_loader.h"
#include "loader/qcovariate_loader.h"
#include "types/fixed_effects.h"

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
        .samples_loaded = phenotype_loader_->sample_ids().size(),
        .trait_name = phenotype_loader_->name()};
}

CovarStats DataPipe::load_covariates()
{
    std::vector<std::string> q_names;
    std::vector<std::string> d_names;

    if (!config_.qcovar_path.empty())
    {
        qcovar_loader_ = std::make_unique<detail::QuantitativeCovariateLoader>(
            config_.qcovar_path, config_.iid_only);
        q_names = qcovar_loader_->column_names();
    }

    if (!config_.dcovar_path.empty())
    {
        dcovar_loader_ = std::make_unique<detail::DiscreteCovariateLoader>(
            config_.dcovar_path, config_.iid_only);
        d_names = dcovar_loader_->column_names();
    }

    return CovarStats{
        .qcovar_loaded = q_names.size(),
        .dcovar_loaded = d_names.size(),
        .q_names = std::move(q_names),
        .d_names = std::move(d_names)};
}

GrmStats DataPipe::load_additive_grm()
{
    if (config_.additive_grm_path.empty())
    {
        throw ArgumentValidationException("Additive GRM path is required.");
    }

    additive_grm_loader_
        = std::make_unique<detail::GrmLoader>(config_.additive_grm_path);

    return GrmStats{
        .samples_in_file
        = static_cast<size_t>(additive_grm_loader_->num_samples())};
}

GrmStats DataPipe::load_dominance_grm()
{
    if (config_.dominance_grm_path.empty())
    {
        throw ArgumentValidationException("Dominance GRM path is required.");
    }

    dominance_grm_loader_
        = std::make_unique<detail::GrmLoader>(config_.dominance_grm_path);

    return GrmStats{
        .samples_in_file
        = static_cast<size_t>(dominance_grm_loader_->num_samples())};
}

IntersectionStats DataPipe::intersect_samples()
{
    size_t total_before = sample_manager_->num_common_samples();

    if (phenotype_loader_)
    {
        total_before
            = std::max(total_before, phenotype_loader_->sample_ids().size());
        sample_manager_->intersect(phenotype_loader_->sample_ids());
    }

    if (qcovar_loader_)
    {
        total_before
            = std::max(total_before, qcovar_loader_->sample_ids().size());
        sample_manager_->intersect(qcovar_loader_->sample_ids());
    }

    if (dcovar_loader_)
    {
        total_before
            = std::max(total_before, dcovar_loader_->sample_ids().size());
        sample_manager_->intersect(dcovar_loader_->sample_ids());
    }

    if (additive_grm_loader_)
    {
        total_before = std::max(
            total_before,
            static_cast<size_t>(additive_grm_loader_->num_samples()));
        sample_manager_->intersect(additive_grm_loader_->sample_ids());
    }

    if (dominance_grm_loader_)
    {
        total_before = std::max(
            total_before,
            static_cast<size_t>(dominance_grm_loader_->num_samples()));
        sample_manager_->intersect(dominance_grm_loader_->sample_ids());
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
    const auto& id_map = sample_manager_->common_id_map();

    phenotype_ = phenotype_loader_ ? phenotype_loader_->load(id_map)
                                   : Eigen::VectorXd::Zero();

    std::optional<QuantitativeCovariate> qcov;
    std::optional<DiscreteCovariate> dcov;

    if (qcovar_loader_)
    {
        qcov = qcovar_loader_->load(id_map);
    }

    if (dcovar_loader_)
    {
        dcov = dcovar_loader_->load(id_map);
    }
    if (!dcov && !qcov)
    {
        fixed_effects_ = FixedEffect::build_bayes(phenotype_.size());
    }
    else
    {
        fixed_effects_ = FixedEffect::build(std::move(qcov), std::move(dcov));
    }

    if (additive_grm_loader_)
    {
        additive_grm_ = std::make_unique<Eigen::MatrixXd>(
            additive_grm_loader_->load(id_map));
    }

    if (dominance_grm_loader_)
    {
        dominance_grm_ = std::make_unique<Eigen::MatrixXd>(
            dominance_grm_loader_->load(id_map));
    }
}

const std::vector<std::string>& DataPipe::fixed_effect_names() const
{
    return fixed_effects_.names;
}

}  // namespace gelex
