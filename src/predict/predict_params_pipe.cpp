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

#include "predict_params_pipe.h"

#include <filesystem>

#include "gelex/exception.h"

namespace gelex
{

PredictParamsPipe::PredictParamsPipe(const Config& config)
{
    if (config.snp_effect_path.empty())
    {
        throw InvalidInputException("SNP effect path must be provided");
    }
    load_snp_effects(config.snp_effect_path);

    if (config.covar_effect_path.empty())
    {
        throw InvalidInputException("params effect path must be provided");
    }
    load_covar_effects(config.covar_effect_path);
}

void PredictParamsPipe::load_snp_effects(const std::filesystem::path& path)
{
    detail::SnpEffectLoader loader(path);
    snp_effects_ = std::move(loader).take_effects();
}

void PredictParamsPipe::load_covar_effects(const std::filesystem::path& path)
{
    detail::CovarEffectLoader loader(path);
    covar_effects_ = std::move(loader).take_effects();
}

}  // namespace gelex
