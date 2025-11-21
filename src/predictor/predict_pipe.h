#pragma once

#include <expected>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "../src/data/loader.h"
#include "Eigen/Core"
#include "covariate_loader.h"
#include "gelex/error.h"

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
        std::filesystem::path qcovar_path;
        std::filesystem::path covar_path;
        bool iid_only = false;
        std::string output_prefix;
    };

    static auto create(
        const Config& config,
        std::shared_ptr<SampleManager> sample_manager)
        -> std::expected<PredictDataPipe, Error>;

    const Eigen::MatrixXd& qcovariates() const { return qcovariates_; }
    Eigen::MatrixXd&& take_qcovariates() && { return std::move(qcovariates_); }

    const std::vector<std::string>& qcovariate_names() const;
    const std::vector<std::string>& covariate_names() const;

    size_t num_qcovariates() const { return qcovariate_names_.size(); }
    size_t num_covariates() const { return covariate_names_.size(); }

   private:
    PredictDataPipe() = default;

    auto load_qcovariates(const Config& config) -> std::expected<void, Error>;
    auto load_covariates(const Config& config) -> std::expected<void, Error>;

    void intersect();
    void format_covariates();

    std::unique_ptr<detail::QcovarLoader> qcovar_loader_;
    std::unique_ptr<detail::CovarPredictLoader> covar_loader_;

    Eigen::MatrixXd qcovariates_;

    std::shared_ptr<SampleManager> sample_manager_;

    std::vector<std::string> qcovariate_names_;
    std::vector<std::string> covariate_names_;
};

}  // namespace gelex
