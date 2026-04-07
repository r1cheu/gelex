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

#include "gelex/pipeline/simulation_engine.h"

#include <optional>
#include <random>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/sim/effect_sampler.h"
#include "gelex/algo/sim/genetic_value_calculator.h"
#include "gelex/algo/sim/phenotype_generator.h"
#include "gelex/data/reader/bim_reader.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/simulation_writer.h"

namespace gelex
{
namespace
{

struct DominanceSpec
{
    double d2;
    bool has_dominance;
    std::vector<EffectSizeClass> effect_classes;
    std::optional<double> positive_prob;
};

auto make_dominance_spec(const SimulationEngine::Config& config)
    -> DominanceSpec
{
    const double d2 = config.dom_heritability.value_or(0.0);
    const bool has_dominance = d2 > 0.0;

    return {
        .d2 = d2,
        .has_dominance = has_dominance,
        .effect_classes = has_dominance ? config.dom_effect_classes
                                        : std::vector<EffectSizeClass>{},
        .positive_prob
        = has_dominance ? config.dom_positive_prob : std::nullopt,
    };
}

}  // namespace

SimulationEngine::SimulationEngine(Config config) : config_(std::move(config))
{
}

auto SimulationEngine::resolve_output_path(
    const std::filesystem::path& output_path,
    const std::filesystem::path& bed_path) -> std::filesystem::path
{
    return output_path.empty() ? std::filesystem::path(bed_path) : output_path;
}

auto SimulationEngine::run(const SimulateObserver& observer) -> void
{
    std::mt19937_64 rng(config_.seed);

    auto bim_path = config_.bed_path;
    bim_path.replace_extension(".bim");
    detail::BimReader bim_reader(bim_path);
    auto snp_ids = bim_reader.get_ids();

    const auto dominance = make_dominance_spec(config_);

    EffectSampler effect_sampler(
        config_.add_effect_classes,
        dominance.effect_classes,
        rng,
        dominance.positive_prob);

    auto causal_effects
        = effect_sampler.sample(static_cast<Eigen::Index>(snp_ids.size()));

    GeneticValueCalculator calculator(
        config_.bed_path, dominance.has_dominance);
    auto genetic_values = calculator.calculate(causal_effects, observer);

    // Generate phenotypes
    PhenotypeGenerator phenotype_generator(
        config_.add_heritability, dominance.d2, config_.intercept, rng);

    auto result = phenotype_generator.generate(genetic_values);

    std::optional<double> actual_dom_positive_prob;
    if (dominance.has_dominance && dominance.positive_prob)
    {
        const auto& dom = causal_effects.dominance;
        Eigen::Index n_positive
            = (dom.array() > 0.0).cast<Eigen::Index>().sum();
        Eigen::Index n_nonzero
            = (dom.array() != 0.0).cast<Eigen::Index>().sum();
        if (n_nonzero > 0)
        {
            actual_dom_positive_prob = static_cast<double>(n_positive)
                                       / static_cast<double>(n_nonzero);
        }
    }

    notify(
        observer,
        HeritabilityGeneratedEvent{
            .additive = result.true_h2,
            .dominance = result.true_d2,
            .dom_positive_prob = actual_dom_positive_prob,
        });

    // Update causal effects with the dominance scaling factor
    auto dom_scale = result.dom_scale;
    if (dom_scale != 1.0)
    {
        causal_effects.dominance *= dom_scale;
    }

    auto output_prefix
        = resolve_output_path(config_.output_path, config_.bed_path);
    SimulationWriter writer(output_prefix);

    writer.write_phenotypes(result.phenotypes, calculator.sample_ids());
    writer.write_causal_effects(snp_ids, causal_effects);

    notify(
        observer,
        OutputsWrittenEvent{
            .phenotype_path = writer.phenotype_path().string(),
            .snp_effect_path = writer.causal_path().string(),
        });
}
}  // namespace gelex
