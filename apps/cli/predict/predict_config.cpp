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

#include <argparse.h>

#include "gelex/data/genotype/bed_path.h"

namespace gelex::cli
{

auto make_predict_config(argparse::ArgumentParser& cmd) -> PredictEngine::Config
{
    return PredictEngine::Config{
        .bed_path = gelex::format_bed_path(cmd.get("bfile")),
        .snp_effect_path = cmd.get("--snp-eff"),
        .covar_effect_path = cmd.get("--covar-eff"),
        .qcovar_path = cmd.get("--qcovar"),
        .dcovar_path = cmd.get("--dcovar"),
        .output_path = cmd.get("--out")};
}

}  // namespace gelex::cli
