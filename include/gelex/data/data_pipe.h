#pragma once

#include <filesystem>
#include <memory>
#include <string>
#include <variant>
#include <vector>

#include <Eigen/Dense>

#include "Eigen/Core"

#include "gelex/data/genotype_loader.h"
#include "gelex/data/genotype_matrix.h"
#include "gelex/data/genotype_mmap.h"
#include "gelex/data/genotype_pipe.h"
#include "gelex/data/sample_manager.h"

namespace gelex
{
namespace detail
{
class PhenotypeLoader;
class QcovarLoader;
class CCovarLoader;
}  // namespace detail
class DataPipe
{
   public:
    struct Config
    {
        std::filesystem::path phenotype_path;
        int phenotype_column = 3;
        std::filesystem::path bed_path;
        bool use_dominance_effect = false;
        bool use_mmap = false;
        int chunk_size = 10000;
        std::filesystem::path qcovar_path;
        std::filesystem::path covar_path;
        bool iid_only = false;
        std::string output_prefix;
    };

    explicit DataPipe(const Config& config);
    DataPipe(const DataPipe&) = delete;
    DataPipe(DataPipe&&) = default;
    DataPipe& operator=(const DataPipe&) = delete;
    DataPipe& operator=(DataPipe&&) = default;
    ~DataPipe();

    Eigen::VectorXd take_phenotype() && { return std::move(phenotype_); }
    Eigen::MatrixXd take_fixed_effects() &&
    {
        return std::move(fixed_effects_);
    }
    std::variant<GenotypeMap, GenotypeMatrix> take_additive_matrix() &&
    {
        return std::move(*additive_matrix_);
    }
    std::variant<GenotypeMap, GenotypeMatrix> take_dominance_matrix() &&
    {
        return std::move(*dominance_matrix_);
    }
    bool has_dominance_matrix() const { return dominance_matrix_ != nullptr; }

    const std::string& phenotype_name() const;
    const std::vector<std::string>& qcovariate_names() const;
    const std::vector<std::string>& covariate_names() const;
    const std::vector<std::string>& fixed_effect_names() const;

    size_t num_qcovariates() const { return qcovariate_names_.size(); }
    size_t num_covariates() const { return covariate_names_.size(); }
    size_t num_fixed_effects() const { return fixed_effect_names_.size(); }

   private:
    DataPipe() = default;

    void load_phenotype(const Config& config);
    void load_qcovariates(const Config& config);
    void load_covariates(const Config& config);
    void load_additive(const Config& config);
    void load_dominance(const Config& config);

    template <typename Processor, typename TargetPtr>
    void load_genotype_impl(
        const Config& config,
        const std::string& suffix,
        TargetPtr& target)
    {
        auto assign_to_target = [&](auto&& result_value)
        {
            using VariantType = typename TargetPtr::element_type;
            target = std::make_unique<VariantType>(
                std::forward<decltype(result_value)>(result_value));
        };

        if (config.use_mmap)
        {
            std::string file_path = config.output_prefix + suffix;
            auto pipe = gelex::GenotypePipe(
                config.bed_path, sample_manager_, file_path);
            auto result = pipe.template process<Processor>(config.chunk_size);
            assign_to_target(std::move(result));
        }
        else
        {
            auto loader
                = gelex::GenotypeLoader(config.bed_path, sample_manager_);
            auto result = loader.template process<Processor>(config.chunk_size);
            assign_to_target(std::move(result));
        }
    }

    void intersect();
    void convert_to_matrices();

    std::unique_ptr<detail::PhenotypeLoader> phenotype_loader_;
    std::unique_ptr<detail::QcovarLoader> qcovar_loader_;
    std::unique_ptr<detail::CCovarLoader> covar_loader_;

    Eigen::VectorXd phenotype_;
    Eigen::MatrixXd fixed_effects_;

    std::shared_ptr<SampleManager> sample_manager_;

    std::unique_ptr<std::variant<GenotypeMap, GenotypeMatrix>> additive_matrix_;
    std::unique_ptr<std::variant<GenotypeMap, GenotypeMatrix>>
        dominance_matrix_;

    std::string phenotype_name_;
    std::vector<std::string> qcovariate_names_;
    std::vector<std::string> covariate_names_;
    std::vector<std::string> fixed_effect_names_;
};

}  // namespace gelex
