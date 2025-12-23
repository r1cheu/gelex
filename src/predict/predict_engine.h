#ifndef GELEX_PREDICT_PREDICT_ENGINE_H
#define GELEX_PREDICT_PREDICT_ENGINE_H

#include <filesystem>
#include <memory>
#include <vector>

#include <Eigen/Core>
#include "../src/data/loader/snp_effect_loader.h"
#include "../src/predict/covar_effect_loader.h"

namespace gelex
{

class PredictParamsPipe;
class PredictDataPipe;

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
        bool iid_only = false;

        void validate() const;
    };

    explicit PredictEngine(const Config& config);

    void run();

    const Eigen::VectorXd& predictions() const { return predictions_; }
    const std::vector<std::string>& sample_ids() const { return sample_ids_; }

   private:
    void load_parameters();
    void load_data();
    void validate_dimensions();
    void compute_predictions();
    void write_output();

    Config config_;
    Eigen::VectorXd predictions_;
    std::vector<std::string> sample_ids_;

    std::unique_ptr<PredictParamsPipe> params_pipe_;
    std::unique_ptr<PredictDataPipe> data_pipe_;

    SnpEffects snp_effects_;
    detail::CovarEffects covar_effects_;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_PREDICT_ENGINE_H
