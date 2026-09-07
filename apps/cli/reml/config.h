// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_REML_CONFIG_H_
#define APPS_CLI_REML_CONFIG_H_

#include <algorithm>
#include <string>
#include <thread>

#include "cli/common_data.h"
#include "cli/reml_data.h"

namespace cli
{

struct RemlConfig
{
    BaseDataConfig base_data;
    RemlDataConfig random;
    std::string out_prefix{"gelex"};
    int max_iter{100};
    double tolerance{1e-6};
    int threads{
        std::max(1, static_cast<int>(std::thread::hardware_concurrency() / 2))};
};

}  // namespace cli

#endif  // APPS_CLI_REML_CONFIG_H_
