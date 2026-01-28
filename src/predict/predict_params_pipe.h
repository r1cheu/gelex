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

#ifndef GELEX_PREDICT_PREDICT_PARAMS_PIPE_H
#define GELEX_PREDICT_PREDICT_PARAMS_PIPE_H

#include <filesystem>

#include "../data/loader/snp_effect_loader.h"
#include "covar_effect_loader.h"

namespace gelex
{

class PredictParamsPipe
{
   public:
    struct Config
    {
        std::filesystem::path snp_effect_path;
        std::filesystem::path covar_effect_path;
    };

    explicit PredictParamsPipe(const Config& config);

    const SnpEffects& snp_effects() const { return snp_effects_; }
    const detail::CovarEffects& covar_effects() const { return covar_effects_; }
    SnpEffects take_snp_effects() && { return std::move(snp_effects_); }
    detail::CovarEffects&& take_covar_effects() &&
    {
        return std::move(covar_effects_);
    }

   private:
    void load_snp_effects(const std::filesystem::path& path);
    void load_covar_effects(const std::filesystem::path& path);

    SnpEffects snp_effects_;
    detail::CovarEffects covar_effects_;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_PREDICT_PARAMS_PIPE_H
