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

#include "gelex/data/genetic_value_calculator.h"

#include <Eigen/Core>

#include "gelex/data/genotype_processor.h"

namespace gelex
{

auto GeneticValueCalculator::calculate_chunk(
    const Eigen::Ref<const Eigen::MatrixXd>& chunk,
    const std::unordered_map<Eigen::Index, CausalEffect>& effects,
    Eigen::Index chunk_start,
    Eigen::Index chunk_end,
    bool has_dominance) -> GeneticValues
{
    const auto n_individuals = chunk.rows();

    Eigen::VectorXd additive_values = Eigen::VectorXd::Zero(n_individuals);
    Eigen::VectorXd dominance_values = Eigen::VectorXd::Zero(n_individuals);

    Eigen::MatrixXd add_chunk = chunk;
    process_matrix<grm::Standardized::Additive>(add_chunk);

    Eigen::MatrixXd dom_chunk;
    if (has_dominance)
    {
        dom_chunk = chunk;
        process_matrix<grm::Standardized::Dominant>(dom_chunk);
    }

    for (const auto& [global_idx, effect] : effects)
    {
        if (global_idx < chunk_start || global_idx >= chunk_end)
        {
            continue;
        }
        const Eigen::Index local_idx = global_idx - chunk_start;

        additive_values += add_chunk.col(local_idx) * effect.additive;

        if (has_dominance)
        {
            dominance_values += dom_chunk.col(local_idx) * effect.dominance;
        }
    }

    return {
        .additive = std::move(additive_values),
        .dominance = std::move(dominance_values)};
}

}  // namespace gelex
