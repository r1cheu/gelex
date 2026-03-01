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

#ifndef GELEX_ALGO_SIM_GENETIC_VALUE_CALCULATOR_H_
#define GELEX_ALGO_SIM_GENETIC_VALUE_CALCULATOR_H_

#include <filesystem>
#include <memory>
#include <span>
#include <string>

#include <Eigen/Core>

#include "gelex/algo/sim/effect_sampler.h"
#include "gelex/data/genotype/bed_pipe.h"
#include "gelex/data/genotype/sample_manager.h"
#include "gelex/infra/logging/simulate_event.h"

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
    GeneticValueCalculator(
        const std::filesystem::path& bed_path,
        bool has_dominance);

    auto calculate(
        const CausalEffects& effects,
        const SimulateObserver& observer = {}) const -> GeneticValues;

    [[nodiscard]] auto sample_ids() const -> std::span<const std::string>;

   private:
    auto encode_chunk(const Eigen::Ref<const Eigen::MatrixXd>& chunk) const
        -> std::pair<Eigen::MatrixXd, Eigen::MatrixXd>;

    bool has_dominance_;
    std::shared_ptr<SampleManager> sample_manager_;
    BedPipe bed_pipe_;

    static constexpr Eigen::Index SNP_CHUNK_SIZE = 10000;
};

}  // namespace gelex

#endif  // GELEX_ALGO_SIM_GENETIC_VALUE_CALCULATOR_H_
