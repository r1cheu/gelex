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

#ifndef GELEX_DATA_SIMULATE_H_
#define GELEX_DATA_SIMULATE_H_

#include <filesystem>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/effect_sampler.h"

namespace gelex
{

class PhenotypeSimulator
{
   public:
    struct Config
    {
        std::filesystem::path bed_path;
        double add_heritability;
        double dom_heritability;
        std::vector<EffectSizeClass> add_effect_classes;
        std::vector<EffectSizeClass> dom_effect_classes;
        double intercept;
        int seed;
        std::filesystem::path output_path;
    };

    explicit PhenotypeSimulator(Config config);

    void simulate();

   private:
    static auto resolve_output_path(
        const std::filesystem::path& output_path,
        const std::filesystem::path& bed_path) -> std::filesystem::path;

    Config config_;

    static constexpr Eigen::Index SNP_CHUNK_SIZE = 10000;
};

}  // namespace gelex

#endif  // GELEX_DATA_SIMULATE_H_
