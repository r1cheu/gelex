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

#ifndef GELEX_ENGINE_GRM_H_
#define GELEX_ENGINE_GRM_H_

#include <string>
#include <vector>

#include "gelex/data/genotype_method.h"
#include "gelex/infra/logging/grm_event.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class GrmEngine
{
   public:
    struct Config
    {
        std::string bfile_prefix;

        std::vector<GeneticMode> requested_effects;
        GenotypeMethod method;
        bool do_loco;

        std::string out_prefix;
        int chunk_size;
    };

    explicit GrmEngine(Config config);

    auto compute(const GrmObserver& observer = {}) -> void;

   private:
    Config config_;
};

}  // namespace gelex

#endif  // GELEX_ENGINE_GRM_H_
