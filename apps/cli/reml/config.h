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

#ifndef APPS_CLI_REML_CONFIG_H_
#define APPS_CLI_REML_CONFIG_H_

#include <algorithm>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include "cli/common_data.h"

namespace cli
{

struct RemlConfig
{
    BaseDataConfig base_data;
    std::vector<std::string> grm_prefixes;
    std::optional<std::string> rand_path;
    std::string out_prefix{"gelex"};
    bool loco{false};
    int max_iter{100};
    double tolerance{1e-6};
    int threads{
        std::max(1, static_cast<int>(std::thread::hardware_concurrency() / 2))};
};

}  // namespace cli

#endif  // APPS_CLI_REML_CONFIG_H_
