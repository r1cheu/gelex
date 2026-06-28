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

#ifndef GELEX_PREDICT_COMPUTE_H_
#define GELEX_PREDICT_COMPUTE_H_

#include <Eigen/Core>

#include "gelex/io/predict/input_reader.h"
#include "gelex/predict/types.h"

namespace gelex::predict
{

[[nodiscard]] auto compute_gebv(
    const GenotypeData& geno,
    const SnpEffects& effects) -> GEBVResult;

[[nodiscard]] auto compute_covariate_effects(
    const Eigen::MatrixXd& covariates,
    const Coefficients& coefficients) -> CovariateResult;

}  // namespace gelex::predict

#endif  // GELEX_PREDICT_COMPUTE_H_
