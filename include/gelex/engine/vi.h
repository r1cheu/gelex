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

#ifndef GELEX_ENGINE_VI_H_
#define GELEX_ENGINE_VI_H_

#include <string>
#include <vector>

#include "gelex/algo/infer/params.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/model/bayes/legacy_method.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel;

namespace vi
{

class Engine
{
   public:
    struct Config
    {
        std::string bfile_prefix;
        bayes::LegacyBayesConfig method;
        std::vector<GeneticMode> requested_effects;
        vi::Params params;
        std::string out_prefix;
    };

    explicit Engine(Config config);
    auto run(
        const BayesModel& model,
        bayes::LegacyBayesMethod method,
        const VIObserver& observer = {}) -> void;

   private:
    Config config_;
};

}  // namespace vi

}  // namespace gelex

#endif  // GELEX_ENGINE_VI_H_
