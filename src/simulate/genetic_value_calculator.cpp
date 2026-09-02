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

#include <Eigen/Core>
#include <cstddef>
#include <filesystem>
#include <ranges>
#include <span>
#include <string>
#include <vector>

#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/notify.h"
#include "gelex/simulate/progress.h"
#include "gelex/simulate/sim_types.h"

namespace gelex
{

GeneticValueCalculator::GeneticValueCalculator(
    const std::filesystem::path& bed_path)
    : bed_(open_bed(bed_path.string()))
{
}

template <GeneticMode Mode>
auto GeneticValueCalculator::calculate(
    GeneticValues& genetic_values,
    GenotypeMethod geno_method,
    const SimulateObserver& observer) const -> Eigen::VectorXd
{
    const auto& causal_snps = genetic_values.causal_snps;
    const Eigen::Index n_individuals = bed_.num_samples();
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
              {
                  return static_cast<Eigen::Index>(bed_.snp_index().at(snp.id));
              })
          | std::ranges::to<std::vector>();

    notify(
        observer,
        SimulateProgressEvent{
            .total = static_cast<std::size_t>(n_causal),
            .current = 0,
            .done = false,
        });

    const LocusEncoder encoder{bed_};
    const EncodingSpec spec{encoding_spec_from_method(Mode, geno_method)};
    Eigen::MatrixXd genotype(n_individuals, n_causal);
    for (const auto [i, variant] : std::views::enumerate(col_indices))
    {
        const LocusStats stats{encoder.count(variant)};
        encoder.expand(
            variant, encoder.encoding(variant, stats, spec), genotype.col(i));
    }

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
    GenotypeMethod,
    const SimulateObserver&) const -> Eigen::VectorXd;

template auto GeneticValueCalculator::calculate<GeneticMode::D>(
    GeneticValues&,
    GenotypeMethod,
    const SimulateObserver&) const -> Eigen::VectorXd;

auto GeneticValueCalculator::sample_ids() const -> std::span<const std::string>
{
    return bed_.sample_index().keys();
}

auto GeneticValueCalculator::snp_ids() const -> std::span<const std::string>
{
    return bed_.snp_index().keys();
}

}  // namespace gelex
