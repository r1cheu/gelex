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

#ifndef GELEX_SIMULATE_GENETIC_VALUE_SCALER_H_
#define GELEX_SIMULATE_GENETIC_VALUE_SCALER_H_

#include <Eigen/Core>
#include <optional>

#include "gelex/simulate/sim_types.h"

namespace gelex
{

class GeneticValueScaler
{
   public:
    GeneticValueScaler(std::optional<double> h2, std::optional<double> d2)
        : h2_(h2), d2_(d2)
    {
    }

    auto scale(
        GeneticValues* additive_values,
        GeneticValues* dominance_values,
        Eigen::Ref<Eigen::VectorXd> residual) const -> void;

   private:
    std::optional<double> h2_;
    std::optional<double> d2_;
};

}  // namespace gelex

#endif  // GELEX_SIMULATE_GENETIC_VALUE_SCALER_H_
