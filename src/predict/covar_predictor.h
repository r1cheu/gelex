#ifndef GELEX_PREDICT_COVAR_PREDICTOR_H
#define GELEX_PREDICT_COVAR_PREDICTOR_H

#include <Eigen/Core>
#include <string>
#include <vector>

#include "covar_effect_loader.h"
#include "predict_pipe.h"

namespace gelex
{

class CovarPredictor
{
   public:
    struct Result
    {
        Eigen::MatrixXd predictions;
        std::vector<std::string> names;
    };

    explicit CovarPredictor(const detail::CovarEffects& effects);

    Result compute(const PredictData& data);

   private:
    void validate_intercept() const;
    void validate_qcovariates() const;
    void validate_continuous_coefficients() const;
    void validate_categorical_coefficients() const;

    void compute_intercept();
    void compute_continuous();
    void compute_categorical();

    const detail::CovarEffects* effects_;
    const PredictData* data_;
    Eigen::Index n_samples_;
    size_t n_cont_;
    size_t n_cat_;
    Result result_;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_COVAR_PREDICTOR_H
