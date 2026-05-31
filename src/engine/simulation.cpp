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

#include "gelex/engine/simulation.h"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/data/reader.h"
#include "gelex/data/sample_id.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/logging/simulate_event.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/io/detail/text_writer.h"
#include "gelex/simulate/effect_sampler.h"
#include "gelex/simulate/genetic_value_calculator.h"
#include "gelex/simulate/genetic_value_scaler.h"
#include "gelex/simulate/sim_types.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

SimulationEngine::SimulationEngine(Config config) : config_(std::move(config))
{
}

auto SimulationEngine::run(const SimulateObserver& observer) -> void
{
    std::mt19937_64 rng(config_.seed);

    auto bim = read_bim(config_.bfile_prefix + ".bim");
    auto fam = read_fam(config_.bfile_prefix + ".fam");

    GeneticValueCalculator calculator(config_.bfile_prefix, bim, fam);

    auto all_ids = bim.index().keys();
    std::vector<std::string_view> shuffled_ids(all_ids.begin(), all_ids.end());
    std::ranges::shuffle(shuffled_ids, rng);

    std::optional<GeneticValues> additive;
    if (config_.additive)
    {
        additive.emplace();
        additive->causal_snps
            = NormalSampler(shuffled_ids, config_.additive->effect_sizes)(rng);
        additive->gebv = calculator.calculate<GeneticMode::A>(
            *additive, config_.geno_method, observer);
    }

    std::optional<GeneticValues> dominance;
    if (config_.dominance)
    {
        dominance.emplace();
        dominance->causal_snps
            = NormalSampler(shuffled_ids, config_.dominance->effect_sizes)(rng);
        dominance->gebv = calculator.calculate<GeneticMode::D>(
            *dominance, config_.geno_method, observer);
    }

    const Eigen::Index n_samples
        = additive ? additive->gebv.size() : dominance->gebv.size();
    std::normal_distribution<double> normal01(0.0, 1.0);
    Eigen::VectorXd residual = Eigen::VectorXd::NullaryExpr(
        n_samples, [&] { return normal01(rng); });

    auto h2_of = [](const auto& scheme) -> std::optional<double>
    {
        return scheme ? std::optional<double>(scheme->heritability)
                      : std::nullopt;
    };
    GeneticValueScaler scaler(
        h2_of(config_.additive), h2_of(config_.dominance));
    scaler.scale(
        additive ? &*additive : nullptr,
        dominance ? &*dominance : nullptr,
        residual);

    Eigen::VectorXd phenotypes = additive ? additive->gebv : dominance->gebv;
    if (additive && dominance)
    {
        phenotypes += dominance->gebv;
    }
    phenotypes += residual;

    const double var_phen = stats::detail::vecvar(phenotypes);
    SimulateVarianceSummaryEvent summary;
    if (additive && var_phen > 0.0)
    {
        summary.realized_h2 = stats::detail::vecvar(additive->gebv) / var_phen;
    }
    if (dominance && var_phen > 0.0)
    {
        summary.realized_d2 = stats::detail::vecvar(dominance->gebv) / var_phen;
    }
    notify(observer, summary);

    io::detail::TextWriter writer(config_.output_prefix + ".phen");
    writer.write_header({"FID", "IID", "Phenotype"});
    for (Eigen::Index i = 0; i < n_samples; ++i)
    {
        auto [fid, iid] = split_sample_id(fam.index().keys()[i]);
        writer.write(fmt::format("{}\t{}\t{}", fid, iid, phenotypes(i)));
    }

    io::detail::TextWriter effect_writer(config_.output_prefix + ".causal");
    auto write_single
        = [&](std::string_view column, const std::vector<CausalSnp>& snps)
    {
        effect_writer.write_header({"id", std::string(column)});
        for (const auto& snp : snps)
        {
            effect_writer.write(fmt::format("{}\t{}", snp.id, snp.effect));
        }
    };
    if (additive && dominance)
    {
        effect_writer.write_header({"id", "additive", "dominance"});
        const auto& add_snps = additive->causal_snps;
        const auto& dom_snps = dominance->causal_snps;
        for (std::size_t i = 0; i < add_snps.size(); ++i)
        {
            effect_writer.write(
                fmt::format(
                    "{}\t{}\t{}",
                    add_snps[i].id,
                    add_snps[i].effect,
                    dom_snps[i].effect));
        }
    }
    else if (additive)
    {
        write_single("additive", additive->causal_snps);
    }
    else
    {
        write_single("dominance", dominance->causal_snps);
    }
}

}  // namespace gelex
