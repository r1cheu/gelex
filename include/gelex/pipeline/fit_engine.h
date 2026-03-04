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

#ifndef GELEX_PIPELINE_FIT_ENGINE_H_
#define GELEX_PIPELINE_FIT_ENGINE_H_

#include <optional>
#include <string>
#include <vector>

#include "gelex/algo/infer/params.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/types/effects.h"

namespace gelex
{
class DataPipe;
class FitEngine
{
   public:
    struct Config
    {
        std::string bfile_prefix;
        BayesAlphabet method;

        int seed;
        MCMCParams mcmc_params;

        std::optional<std::vector<double>> pi;
        std::optional<std::vector<double>> dpi;
        std::optional<std::vector<double>> scale;
        std::optional<std::vector<double>> dscale;

        std::string out_prefix;
    };

    explicit FitEngine(Config config);
    auto run(DataPipe&& data, const FitObserver& observer = {}) -> void;

   private:
    Config config_;
};

}  // namespace gelex

#endif  // GELEX_PIPELINE_FIT_ENGINE_H_
