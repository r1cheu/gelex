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

#ifndef APPS_CLI_POST_CONFIG_H_
#define APPS_CLI_POST_CONFIG_H_

#include <optional>
#include <string>
#include <vector>

namespace CLI
{
class App;
}

namespace cli
{

struct PostConfig
{
    std::vector<std::string> in_prefixes;
    std::optional<std::string> gfile;
    double hdpi_width;
};

auto make_post_config(CLI::App& cmd) -> PostConfig;

}  // namespace cli

#endif  // APPS_CLI_POST_CONFIG_H_
