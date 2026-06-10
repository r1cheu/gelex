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

#ifndef GELEX_SIMULATE_GENETIC_VALUE_CALCULATOR_H_
#define GELEX_SIMULATE_GENETIC_VALUE_CALCULATOR_H_

#include <filesystem>
#include <span>
#include <string>

#include <Eigen/Core>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype/bed_pipe.h"
#include "gelex/data/genotype/method.h"
#include "gelex/infra/logging/simulate_event.h"
#include "gelex/simulate/sim_types.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class GeneticValueCalculator
{
   public:
    GeneticValueCalculator(
        const std::filesystem::path& bed_path,
        const dataframe::DataFrame<std::string>& bim,
        const dataframe::DataFrame<std::string>& fam);

    template <GeneticMode Mode>
    auto calculate(
        GeneticValues& genetic_values,
        GenotypeMethod geno_method,
        const SimulateObserver& observer = {}) const -> Eigen::VectorXd;

    [[nodiscard]] auto sample_ids() const -> std::span<const std::string>;

   private:
    const dataframe::Index<std::string>* sample_index_;
    const dataframe::Index<std::string>* snp_index_;
    genotype::BedPipe bed_pipe_;
};

}  // namespace gelex

#endif  // GELEX_SIMULATE_GENETIC_VALUE_CALCULATOR_H_
