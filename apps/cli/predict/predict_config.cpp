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

auto make_predict_configs(argparse::ArgumentParser& cmd) -> PredictConfigs
{
    using std::filesystem::path;

    PredictConfigs configs;

    auto gfile = cmd.get<std::string>("--gfile");

    configs.params.snp_effect_path = gfile + ".snp.eff";
    configs.params.covar_effect_path = gfile + ".covar.eff";
    configs.params.sbin_path = gfile + ".sbin";

    configs.data.bed_path = gelex::format_bed_path(cmd.get("bfile"));
    configs.data.qcovar_path
        = cmd.is_used("--qcovar")
              ? std::make_optional<path>(cmd.get("--qcovar"))
              : std::nullopt;
    configs.data.dcovar_path
        = cmd.is_used("--dcovar")
              ? std::make_optional<path>(cmd.get("--dcovar"))
              : std::nullopt;

    configs.engine.output_path = cmd.get("--out");
    configs.engine.chunk_size = cmd.get<int>("--chunk-size");

    return configs;
}

}  // namespace gelex::cli
