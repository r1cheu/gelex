// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_PREDICT_IO_H_
#define APPS_CLI_PREDICT_IO_H_

#include <Eigen/Core>
#include <filesystem>
#include <span>
#include <string>

#include "gelex/genetic_mode.h"

#include "compute.h"

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
