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

#ifndef GELEX_PREDICT_TYPES_H_
#define GELEX_PREDICT_TYPES_H_

#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/io/locistats/reader.h"

namespace gelex::predict
{

struct SbinData
{
    LociStats add;
    LociStats dom;
    bool has_dom{false};
};

struct GenotypeData
{
    Eigen::MatrixXd add;
    std::optional<Eigen::MatrixXd> dom;
};

struct SnpEffects
{
    Eigen::VectorXd add;
    std::optional<Eigen::VectorXd> dom;
};

struct GEBVResult
{
    Eigen::VectorXd total;
    Eigen::VectorXd add_predictions;
    std::optional<Eigen::VectorXd> dom_predictions;
};

struct CovariateResult
{
    Eigen::VectorXd total;
    Eigen::MatrixXd per_covariate;
    std::vector<std::string> covar_names;
};

struct PredictResult
{
    std::vector<std::string> sample_ids;
    Eigen::VectorXd predictions;
    Eigen::VectorXd snp_predictions;
    Eigen::VectorXd add_predictions;
    std::optional<Eigen::VectorXd> dom_predictions;
    Eigen::MatrixXd covar_predictions;
    std::vector<std::string> covar_names;
};

}  // namespace gelex::predict

#endif  // GELEX_PREDICT_TYPES_H_
