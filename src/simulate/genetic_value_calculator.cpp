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

#include "gelex/simulate/genetic_value_calculator.h"

#include <cstddef>
#include <filesystem>
#include <ranges>
#include <span>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/genotype/process_method.h"
#include "gelex/data/genotype/processor.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/logging/simulate_event.h"
#include "gelex/simulate/sim_types.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

GeneticValueCalculator::GeneticValueCalculator(
    const std::filesystem::path& bed_path,
    const dataframe::DataFrame<std::string>& bim,
    const dataframe::DataFrame<std::string>& fam)
    : sample_index_(&fam.index()),
      snp_index_(&bim.index()),
      bed_pipe_(bed_path, *sample_index_)
{
}

template <GeneticMode Mode>
auto GeneticValueCalculator::calculate(
    GeneticValues& genetic_values,
    GenotypeProcessMethod geno_method,
    const SimulateObserver& observer) const -> Eigen::VectorXd
{
    const auto& causal_snps = genetic_values.causal_snps;
    const Eigen::Index n_individuals = bed_pipe_.num_samples();
    const auto n_causal = static_cast<Eigen::Index>(causal_snps.size());

    if (n_causal == 0)
    {
        genetic_values.coeff.resize(0);
        return Eigen::VectorXd::Zero(n_individuals);
    }

    auto col_indices
        = causal_snps
          | std::views::transform(
              [this](const auto& snp)
              { return static_cast<Eigen::Index>(snp_index_->at(snp.id)); })
          | std::ranges::to<std::vector>();

    notify(
        observer,
        SimulateProgressEvent{
            .total = static_cast<std::size_t>(n_causal),
            .current = 0,
            .done = false,
        });

    Eigen::MatrixXd genotype = bed_pipe_.select(col_indices);
    genotype::process_matrix<Mode>(geno_method, genotype);

    genetic_values.coeff.resize(n_causal);
    for (Eigen::Index i = 0; i < n_causal; ++i)
    {
        genetic_values.coeff(i)
            = causal_snps[static_cast<std::size_t>(i)].effect;
    }

    Eigen::VectorXd values = genotype * genetic_values.coeff;

    notify(
        observer,
        SimulateProgressEvent{
            .total = static_cast<std::size_t>(n_causal),
            .current = static_cast<std::size_t>(n_causal),
            .done = true,
        });

    return values;
}

template auto GeneticValueCalculator::calculate<GeneticMode::A>(
    GeneticValues&,
    GenotypeProcessMethod,
    const SimulateObserver&) const -> Eigen::VectorXd;

template auto GeneticValueCalculator::calculate<GeneticMode::D>(
    GeneticValues&,
    GenotypeProcessMethod,
    const SimulateObserver&) const -> Eigen::VectorXd;

auto GeneticValueCalculator::sample_ids() const -> std::span<const std::string>
{
    return sample_index_->keys();
}

}  // namespace gelex
