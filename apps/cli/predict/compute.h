// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
