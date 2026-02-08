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

#ifndef GELEX_CLI_GRM_ARGS_H_
#define GELEX_CLI_GRM_ARGS_H_

#include <filesystem>
#include <string>

namespace argparse
{
class ArgumentParser;
}

struct GrmConfig
{
    std::filesystem::path bed_path;
    std::string out_prefix;
    std::string method;
    int chunk_size;
    bool do_additive;
    bool do_dominant;
    bool do_loco;
    int threads;
};

void setup_grm_args(argparse::ArgumentParser& cmd);

#endif  // GELEX_CLI_GRM_ARGS_H_
