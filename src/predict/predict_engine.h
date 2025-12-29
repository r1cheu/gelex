#ifndef GELEX_PREDICT_PREDICT_ENGINE_H
#define GELEX_PREDICT_PREDICT_ENGINE_H

#include <filesystem>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader/snp_effect_loader.h"
#include "../src/predict/covar_effect_loader.h"
#include "predict_pipe.h"

namespace gelex
{

class PredictWriter;

class PredictEngine
{
   public:
    struct Config
    {
        std::filesystem::path bed_path;
        std::filesystem::path snp_effect_path;
        std::filesystem::path covar_effect_path;
        std::filesystem::path qcovar_path;
        std::filesystem::path dcovar_path;
        std::filesystem::path output_path;
        bool iid_only = false;

        void validate() const;
    };

    explicit PredictEngine(const Config& config);

    void run();

    const Eigen::VectorXd& predictions() const { return predictions_; }
    const std::vector<std::string>& sample_ids() const { return sample_ids_; }
    const Eigen::VectorXd& snp_predictions() const { return snp_predictions_; }
    const Eigen::VectorXd& add_predictions() const { return add_predictions_; }
    const Eigen::VectorXd& dom_predictions() const { return dom_predictions_; }
    const Eigen::MatrixXd& covar_predictions() const
    {
        return covar_predictions_;
    }
    const std::vector<std::string>& covar_prediction_names() const
    {
        return covar_prediction_names_;
    }

   private:
    void load_parameters();
    void load_data();
    void validate_dimensions();
    void compute();
    void write();

    Config config_;
    Eigen::VectorXd predictions_;
    Eigen::VectorXd snp_predictions_;
    Eigen::VectorXd add_predictions_;
    Eigen::VectorXd dom_predictions_;

    std::vector<std::string> sample_ids_;

    Eigen::MatrixXd covar_predictions_;
    std::vector<std::string> covar_prediction_names_;

    PredictData data_;
    SnpEffects snp_effects_;
    detail::CovarEffects covar_effects_;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_PREDICT_ENGINE_H
