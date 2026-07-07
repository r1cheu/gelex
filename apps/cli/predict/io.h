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
#include <span>
#include <string>

#include <Eigen/Core>

#include "compute.h"
#include "gelex/types/genetic_mode.h"

namespace cli
{

auto write_predictions(
    const std::filesystem::path& output_path,
    std::span<const std::string> sample_ids,
    const Eigen::Ref<const Eigen::VectorXd>& prediction,
    const CovariateResult& covar,
    const gelex::ModeMap<Eigen::VectorXd>& gebvs) -> void;

}  // namespace cli

#endif  // APPS_CLI_PREDICT_IO_H_
