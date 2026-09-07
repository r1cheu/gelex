// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_PREDICT_CONFIG_H_
#define APPS_CLI_PREDICT_CONFIG_H_

#include <optional>
#include <string>

namespace cli
{

struct PredictConfig
{
    std::string bfile;
    std::string gfile;
    std::optional<std::string> qcovar;
    std::optional<std::string> dcovar;
    std::string out;
    int chunk_size{10000};
};

}  // namespace cli

#endif  // APPS_CLI_PREDICT_CONFIG_H_
