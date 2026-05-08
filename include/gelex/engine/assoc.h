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

#ifndef GELEX_ENGINE_ASSOC_H_
#define GELEX_ENGINE_ASSOC_H_

#include <string>

#include "gelex/algo/gwas/assoc_type.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype/process_method.h"
#include "gelex/infra/logging/assoc_event.h"
#include "gelex/infra/logging/reml_event.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class PhenoPipe;
class GrmPipe;

class AssocEngine
{
   public:
    struct Config
    {
        GeneticMode mode;

        GenotypeProcessMethod method;
        int chunk_size;

        int max_iter;
        double tol;

        std::string bfile_prefix;
        std::string out_prefix;

        bool loco = false;
        AssocType test_type = AssocType::Single;
    };

    explicit AssocEngine(Config config);

    auto run(
        PhenoPipe& pheno,
        GrmPipe& grm,
        const dataframe::Index<std::string>& sample_index,
        const AssocObserver& observer = {},
        const RemlObserver& reml_observer = {}) -> void;

   private:
    Config config_;
};

}  // namespace gelex

#endif  // GELEX_ENGINE_ASSOC_H_
