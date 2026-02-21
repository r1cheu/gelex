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

#ifndef GELEX_PIPELINE_SIM_PHENOTYPE_SIMULATION_ENGINE_H_
#define GELEX_PIPELINE_SIM_PHENOTYPE_SIMULATION_ENGINE_H_

#include <filesystem>
#include <optional>
#include <vector>

#include "gelex/algo/sim/effect_sampler.h"
#include "gelex/infra/logging/simulate_event.h"

namespace gelex
{

class PhenotypeSimulationEngine
{
   public:
    struct Config
    {
        std::filesystem::path bed_path;
        std::filesystem::path output_path;

        double intercept;

        double add_heritability;
        std::vector<EffectSizeClass> add_effect_classes;

        std::optional<double> dom_heritability;
        std::vector<EffectSizeClass> dom_effect_classes;

        int seed;
    };

    explicit PhenotypeSimulationEngine(Config config);

    auto simulate(const SimulateObserver& observer = {}) -> void;

   private:
    static auto resolve_output_path(
        const std::filesystem::path& output_path,
        const std::filesystem::path& bed_path) -> std::filesystem::path;

    Config config_;
};

}  // namespace gelex

#endif  // GELEX_PIPELINE_SIM_PHENOTYPE_SIMULATION_ENGINE_H_
