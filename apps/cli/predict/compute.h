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

#ifndef APPS_CLI_PREDICT_COMPUTE_H_
#define APPS_CLI_PREDICT_COMPUTE_H_

#include <Eigen/Core>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/encode.h"
#include "gelex/genetic_mode.h"

namespace cli
{

struct CovariateDesign
{
    Eigen::MatrixXd matrix;
    std::vector<std::pair<std::string, gelex::LevelMismatch>> level_mismatches;
};

struct CovariateResult
{
    Eigen::MatrixXd per_covariate;
    std::vector<std::string> covar_names;
};

[[nodiscard]] auto compute_gebv(
    const gelex::ModeMap<Eigen::MatrixXd>& geno,
    const gelex::ModeMap<Eigen::VectorXd>& effects)
    -> gelex::ModeMap<Eigen::VectorXd>;

[[nodiscard]] auto build_covariate_design(
    std::span<const std::string> term_names,
    const std::optional<gelex::DataFrame<std::string>>& qcovar_df,
    const std::optional<gelex::DataFrame<std::string>>& dcovar_df,
    Eigen::Index n_samples) -> CovariateDesign;

[[nodiscard]] auto compute_covariate_effects(
    const Eigen::MatrixXd& covariates,
    std::span<const std::string> term_names,
    const Eigen::Ref<const Eigen::VectorXd>& coefficients) -> CovariateResult;

}  // namespace cli

#endif  // APPS_CLI_PREDICT_COMPUTE_H_
