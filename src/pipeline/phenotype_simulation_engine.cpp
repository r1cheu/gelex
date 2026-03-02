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

#include "gelex/pipeline/phenotype_simulation_engine.h"

#include <random>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/sim/effect_sampler.h"
#include "gelex/algo/sim/genetic_value_calculator.h"
#include "gelex/algo/sim/phenotype_generator.h"
#include "gelex/data/loader/bim_loader.h"
#include "gelex/io/sim/simulation_writer.h"

namespace gelex
{
namespace
{

struct DominanceSpec
{
    double d2;
    bool has_dominance;
    std::vector<EffectSizeClass> effect_classes;
};

auto make_dominance_spec(const PhenotypeSimulationEngine::Config& config)
    -> DominanceSpec
{
    const double d2 = config.dom_heritability.value_or(0.0);
    const bool has_dominance = d2 > 0.0;

    return {
        .d2 = d2,
        .has_dominance = has_dominance,
        .effect_classes = has_dominance ? config.dom_effect_classes
                                        : std::vector<EffectSizeClass>{},
    };
}

}  // namespace

PhenotypeSimulationEngine::PhenotypeSimulationEngine(Config config)
    : config_(std::move(config))
{
}

auto PhenotypeSimulationEngine::resolve_output_path(
    const std::filesystem::path& output_path,
    const std::filesystem::path& bed_path) -> std::filesystem::path
{
    return output_path.empty() ? std::filesystem::path(bed_path) : output_path;
}

auto PhenotypeSimulationEngine::run(const SimulateObserver& observer) -> void
{
    std::mt19937_64 rng(config_.seed);

    auto bim_path = config_.bed_path;
    bim_path.replace_extension(".bim");
    detail::BimLoader bim_loader(bim_path);
    auto snp_ids = bim_loader.get_ids();

    const auto dominance = make_dominance_spec(config_);

    EffectSampler effect_sampler(
        config_.add_effect_classes, dominance.effect_classes, rng);

    auto causal_effects
        = effect_sampler.sample(static_cast<Eigen::Index>(snp_ids.size()));

    GeneticValueCalculator calculator(
        config_.bed_path, dominance.has_dominance);
    auto genetic_values = calculator.calculate(causal_effects, observer);

    // Generate phenotypes
    PhenotypeGenerator phenotype_generator(
        config_.add_heritability, dominance.d2, config_.intercept, rng);

    auto phenotypes = phenotype_generator.generate(genetic_values, observer);

    // Update causal effects with the dominance scaling factor
    double dom_scale = phenotype_generator.dom_scale();
    if (dom_scale != 1.0)
    {
        causal_effects.dominance *= dom_scale;
    }

    auto output_prefix
        = resolve_output_path(config_.output_path, config_.bed_path);
    SimulationWriter writer(output_prefix);

    writer.write_phenotypes(phenotypes, calculator.sample_ids());
    writer.write_causal_effects(snp_ids, causal_effects);

    if (observer)
    {
        SimulateEvent event;
        event.emplace<OutputsWrittenEvent>(OutputsWrittenEvent{
            .phenotype_path = writer.phenotype_path().string(),
            .snp_effect_path = writer.causal_path().string(),
        });
        observer(event);
    }
}
}  // namespace gelex
