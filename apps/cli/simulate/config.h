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

#ifndef APPS_CLI_SIMULATE_CONFIG_H_
#define APPS_CLI_SIMULATE_CONFIG_H_

#include <optional>
#include <string>
#include <vector>

namespace cli
{

struct SimulateConfig
{
    std::string bfile;
    std::string out{"sim.phen"};
    std::optional<double> h2;
    std::vector<double> add_var;
    std::vector<int> add_n;
    std::optional<double> d2;
    std::vector<double> dom_var;
    std::vector<int> dom_n;
    std::optional<double> dom_pos_prob;
    std::string geno_method{"OS"};
    int seed{42};
};

}  // namespace cli

#endif  // APPS_CLI_SIMULATE_CONFIG_H_
