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

#ifndef GELEX_DATA_GENETIC_VALUE_CALCULATOR_H_
#define GELEX_DATA_GENETIC_VALUE_CALCULATOR_H_

#include <unordered_map>

#include <Eigen/Core>

#include "gelex/data/effect_sampler.h"

namespace gelex
{

struct GeneticValues
{
    Eigen::VectorXd additive;
    Eigen::VectorXd dominance;
};

class GeneticValueCalculator
{
   public:
    static auto calculate_chunk(
        const Eigen::Ref<const Eigen::MatrixXd>& chunk,
        const std::unordered_map<Eigen::Index, CausalEffect>& effects,
        Eigen::Index chunk_start,
        Eigen::Index chunk_end,
        bool has_dominance) -> GeneticValues;
};

}  // namespace gelex

#endif  // GELEX_DATA_GENETIC_VALUE_CALCULATOR_H_
