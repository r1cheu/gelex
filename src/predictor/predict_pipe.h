#ifndef GELEX_PREDICTOR_PREDICT_PIPE_H
#define GELEX_PREDICTOR_PREDICT_PIPE_H

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "covariate_loader.h"
#include "data/loader/qcovariate_loader.h"

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

        bool iid_only = false;
    };

    explicit PredictDataPipe(const Config& config);

    auto take_qcovariates() && -> Eigen::MatrixXd
    {
        return std::move(qcovariates_);
    }
    auto take_covariates() && -> std::map<std::string, std::vector<std::string>>
    {
        return std::move(covariates_);
    }
    auto take_genotypes() && -> Eigen::MatrixXd
    {
        return std::move(genotypes_);
    }

    auto qcovariate_names() const -> const std::vector<std::string>&
    {
        return qcovar_loader_->names();
    }
    auto covariate_names() const -> const std::vector<std::string>&
    {
        return covar_loader_->names();
    };

    size_t num_qcovariates() const { return qcovariate_names().size(); }
    size_t num_covariates() const { return covariate_names().size(); }

   private:
    auto load_qcovariates(const Config& config) -> void;
    auto load_covariates(const Config& config) -> void;
    auto load_genotype(const Config& config) -> void;

    void intersect();
    void format_covariates();

    std::unique_ptr<detail::QcovarLoader> qcovar_loader_;
    std::unique_ptr<detail::CovarPredictLoader> covar_loader_;

    Eigen::MatrixXd qcovariates_;
    std::map<std::string, std::vector<std::string>> covariates_;
    Eigen::MatrixXd genotypes_;

    std::shared_ptr<SampleManager> sample_manager_;
};

}  // namespace gelex
#endif  // GELEX_PREDICTOR_PREDICT_PIPE_H
