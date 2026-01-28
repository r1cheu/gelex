/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
