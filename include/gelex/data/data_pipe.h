#pragma once

#include <expected>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "../src/data/loader.h"
#include "Eigen/Core"
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"

namespace gelex
{
class DataPipe
{
   public:
    struct Config
    {
        std::filesystem::path phenotype_path;
        int phenotype_column = 3;
        std::filesystem::path qcovar_path;
        std::filesystem::path covar_path;
        bool iid_only = false;
        std::string output_prefix;
    };

    static auto create(
        const Config& config,
        std::shared_ptr<SampleManager> sample_manager)
        -> std::expected<DataPipe, Error>;

    const Eigen::VectorXd& phenotype() const { return phenotype_; }
    const Eigen::MatrixXd& fixed_effects() const { return fixed_effects_; }

    const std::string& phenotype_name() const;
    const std::vector<std::string>& qcovariate_names() const;
    const std::vector<std::string>& covariate_names() const;
    const std::vector<std::string>& fixed_effect_names() const;

    size_t num_qcovariates() const { return qcovariate_names_.size(); }
    size_t num_covariates() const { return covariate_names_.size(); }
    size_t num_fixed_effects() const { return fixed_effect_names_.size(); }

    bool has_phenotype() const { return !phenotype_.isZero(0); }
    bool has_fixed_effects() const { return fixed_effects_.cols() > 0; }

   private:
    explicit DataPipe(std::shared_ptr<SampleManager> sample_manager)
        : sample_manager_(std::move(sample_manager))
    {
    }

    auto load_phenotype(const Config& config) -> std::expected<void, Error>;
    auto load_qcovariates(const Config& config) -> std::expected<void, Error>;
    auto load_covariates(const Config& config) -> std::expected<void, Error>;

    void intersect();
    void convert_to_matrices();

    std::unique_ptr<detail::PhenotypeLoader> phenotype_loader_;
    std::unique_ptr<detail::QcovarLoader> qcovar_loader_;
    std::unique_ptr<detail::CovarLoader> covar_loader_;

    Eigen::VectorXd phenotype_;
    Eigen::MatrixXd fixed_effects_;

    std::shared_ptr<SampleManager> sample_manager_;

    std::string phenotype_name_;
    std::vector<std::string> qcovariate_names_;
    std::vector<std::string> covariate_names_;
    std::vector<std::string> fixed_effect_names_;
};

}  // namespace gelex
