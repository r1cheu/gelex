#ifndef GELEX_PREDICTOR_PREDICT_PIPE_H
#define GELEX_PREDICTOR_PREDICT_PIPE_H

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader/qcovariate_loader.h"
#include "../src/predictor/predict_dcovariate_loader.h"

namespace gelex
{
class SampleManager;
namespace detail
{
class QcovarLoader;
class DcovarPredictLoader;
}  // namespace detail

class PredictDataPipe
{
   public:
    struct Config
    {
        std::filesystem::path bed_path;
        std::filesystem::path qcovar_path;
        std::filesystem::path dcovar_path;

        bool iid_only = false;
    };

    explicit PredictDataPipe(const Config& config);

    auto take_qcovariates() && -> Eigen::MatrixXd
    {
        return std::move(qcovariates_);
    }
    auto
    take_dcovariates() && -> std::map<std::string, std::vector<std::string>>
    {
        return std::move(dcovariates_);
    }
    auto take_genotypes() && -> Eigen::MatrixXd
    {
        return std::move(genotypes_);
    }

    auto qcovariate_names() const -> const std::vector<std::string>&
    {
        return qcovariate_names_;
    }
    auto dcovariate_names() const -> const std::vector<std::string>&
    {
        return dcovariate_names_;
    };

    size_t num_qcovariates() const { return qcovariate_names_.size(); }
    size_t num_dcovariates() const { return dcovariate_names_.size(); }

   private:
    auto load_qcovariates(const Config& config) -> void;
    auto load_dcovariates(const Config& config) -> void;
    auto load_genotype(const Config& config) -> void;

    void intersect();
    void format_dcovariates();

    std::unique_ptr<detail::QcovarLoader> qcovar_loader_;
    std::unique_ptr<detail::DcovarPredictLoader> dcovar_loader_;
    Eigen::MatrixXd qcovariates_;
    std::vector<std::string> qcovariate_names_;
    std::map<std::string, std::vector<std::string>> dcovariates_;
    std::vector<std::string> dcovariate_names_;

    Eigen::MatrixXd genotypes_;

    std::shared_ptr<SampleManager> sample_manager_;
};

}  // namespace gelex
#endif  // GELEX_PREDICTOR_PREDICT_PIPE_H
