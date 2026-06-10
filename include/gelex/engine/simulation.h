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

#ifndef GELEX_ENGINE_SIMULATION_H_
#define GELEX_ENGINE_SIMULATION_H_

#include <optional>
#include <string>
#include <vector>

#include "gelex/data/genotype/method.h"
#include "gelex/infra/logging/simulate_event.h"
#include "gelex/simulate/sim_types.h"

namespace gelex
{

class SimulationEngine
{
   public:
    struct SimulateScheme
    {
        double heritability;
        std::vector<EffectSize> effect_sizes;
    };

    struct Config
    {
        int seed;
        std::string bfile_prefix;
        std::string output_prefix;

        GenotypeMethod geno_method;

        std::optional<SimulateScheme> additive;
        std::optional<SimulateScheme> dominance;
        std::optional<double> dom_positive_prob;
    };

    explicit SimulationEngine(Config config);

    auto run(const SimulateObserver& observer = {}) -> void;

   private:
    Config config_;
};

}  // namespace gelex

#endif  // GELEX_ENGINE_SIMULATION_H_
