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

#ifndef APPS_CLI_PREDICT_IO_H_
#define APPS_CLI_PREDICT_IO_H_

#include <filesystem>
#include <optional>

#include <Eigen/Core>

#include "gelex/data/dataframe/dataframe.h"
#include "types.h"

namespace cli
{

[[nodiscard]] auto read_coefficients(const std::filesystem::path& path)
    -> Coefficients;

[[nodiscard]] auto read_covariates(
    const std::optional<std::filesystem::path>& qcovar_path,
    const std::optional<std::filesystem::path>& dcovar_path,
    const Coefficients& coefficients,
    gelex::DataFrame<std::string>& sample_df) -> Eigen::MatrixXd;

auto write_predictions(
    const std::filesystem::path& output_path,
    const PredictResult& result) -> void;

}  // namespace cli

#endif  // APPS_CLI_PREDICT_IO_H_
