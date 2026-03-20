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

#include "predict_config.h"

#include <filesystem>
#include <optional>
#include <string>

#include <argparse.h>

#include "gelex/data/genotype/bed_path.h"

namespace gelex::cli
{

auto make_predict_config(argparse::ArgumentParser& cmd) -> PredictEngine::Config
{
    using std::filesystem::path;

    PredictEngine::Config config;

    auto gfile = cmd.get<std::string>("--gfile");
    auto bfile = cmd.get<std::string>("--bfile");

    config.bed_path = gelex::format_bed_path(bfile);
    config.bfile_prefix = bfile;
    config.gfile_prefix = gfile;

    config.qcovar_path = cmd.is_used("--qcovar")
                             ? std::make_optional<path>(cmd.get("--qcovar"))
                             : std::nullopt;
    config.dcovar_path = cmd.is_used("--dcovar")
                             ? std::make_optional<path>(cmd.get("--dcovar"))
                             : std::nullopt;

    config.output_path = cmd.get("--out");

    return config;
}

}  // namespace gelex::cli
