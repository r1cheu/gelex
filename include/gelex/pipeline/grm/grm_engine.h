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

#ifndef GELEX_PIPELINE_GRM_GRM_ENGINE_H_
#define GELEX_PIPELINE_GRM_GRM_ENGINE_H_

#include <filesystem>
#include <string>

#include "gelex/data/genotype/genotype_method_dispatch.h"
#include "gelex/infra/logging/grm_event.h"
#include "gelex/types/freq_effect.h"

namespace gelex
{

class GrmEngine
{
   public:
    struct Config
    {
        std::filesystem::path bed_path;
        std::string out_prefix;
        GenotypeProcessMethod method;
        gelex::freq::GrmType mode;
        int chunk_size;
        bool do_loco;
    };

    explicit GrmEngine(Config config);

    auto compute(const GrmObserver& observer = {}) -> void;

   private:
    Config config_;
};

}  // namespace gelex

#endif  // GELEX_PIPELINE_GRM_GRM_ENGINE_H_
