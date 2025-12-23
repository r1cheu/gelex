#ifndef GELEX_PREDICT_PREDICT_PIPE_H
#define GELEX_PREDICT_PREDICT_PIPE_H

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader/qcovariate_loader.h"
#include "../src/predict/predict_dcovariate_loader.h"

namespace gelex
{
class SampleManager;
namespace detail
{
class QcovarLoader;
class DcovarPredictLoader;
}  // namespace detail

struct PredictData
{
    std::vector<std::string> qcovariate_names;
    Eigen::MatrixXd qcovariates;

    std::vector<std::string> dcovariate_names;
    std::map<std::string, std::vector<std::string>> dcovariates;

    Eigen::MatrixXd genotype;
};

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

    auto take_data() && -> PredictData
    {
        return PredictData{
            .qcovariate_names = std::move(qcovariate_names_),
            .qcovariates = std::move(qcovariates_),
            .dcovariate_names = std::move(dcovariate_names_),
            .dcovariates = std::move(dcovariates_),
            .genotype = std::move(genotypes_)};
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
#endif  // GELEX_PREDICT_PREDICT_PIPE_H
