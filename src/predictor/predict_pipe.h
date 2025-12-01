#pragma once

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

// #include "../src/data/loader.h"  // 不再需要，使用前向声明
#include "../src/data/loader/qcovariate_loader.h"
#include "Eigen/Core"
#include "covariate_loader.h"
#include "snp_effect_loader.h"

namespace gelex
{
class SampleManager;
namespace detail
{
class QcovarLoader;
}

class PredictDataPipe
{
   public:
    struct Config
    {
        std::filesystem::path bed_path;
        std::filesystem::path qcovar_path;
        std::filesystem::path covar_path;

        std::filesystem::path snp_effect_path;
        bool iid_only = false;
        std::string output_prefix;
    };

    explicit PredictDataPipe(const Config& config);

    const Eigen::MatrixXd& qcovariates() const;
    Eigen::MatrixXd take_qcovariates() &&;
    std::map<std::string, std::vector<std::string>> take_covariates() &&;

    const std::vector<std::string>& qcovariate_names() const;
    const std::vector<std::string>& covariate_names() const;

    size_t num_qcovariates() const;
    size_t num_covariates() const;

   private:
    auto load_qcovariates(const Config& config) -> void;
    auto load_covariates(const Config& config) -> void;
    auto load_genotype(const Config& config) -> void;
    auto load_snp_effect(const Config& config) -> void;

    auto process_raw_genotype(const Config& config) -> void;

    void intersect();
    void format_covariates();

    std::unique_ptr<detail::QcovarLoader> qcovar_loader_;
    std::unique_ptr<detail::CovarPredictLoader> covar_loader_;

    Eigen::MatrixXd qcovariates_;
    std::map<std::string, std::vector<std::string>> covariates_;
    SnpEffects snp_effects_;
    std::vector<Eigen::MatrixXd> genotypes_;
    bool has_dom_ = false;

    std::shared_ptr<SampleManager> sample_manager_;

    std::vector<std::string> qcovariate_names_;
    std::vector<std::string> covariate_names_;
};

}  // namespace gelex
