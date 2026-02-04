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

#ifndef GELEX_DATA_SIMULATION_WRITER_H_
#define GELEX_DATA_SIMULATION_WRITER_H_

#include <filesystem>
#include <span>
#include <string>
#include <unordered_map>

#include <Eigen/Core>

#include "gelex/data/effect_sampler.h"

namespace gelex
{

class SimulationWriter
{
   public:
    explicit SimulationWriter(std::filesystem::path output_prefix);

    void write_phenotypes(
        const Eigen::Ref<const Eigen::VectorXd>& phenotypes,
        std::span<const std::string> sample_ids) const;

    void write_causal_effects(
        std::span<const std::string> snp_ids,
        const std::unordered_map<Eigen::Index, CausalEffect>& effects) const;

    [[nodiscard]] auto phenotype_path() const -> std::filesystem::path;
    [[nodiscard]] auto causal_path() const -> std::filesystem::path;

   private:
    std::filesystem::path output_prefix_;
};

}  // namespace gelex

#endif  // GELEX_DATA_SIMULATION_WRITER_H_
