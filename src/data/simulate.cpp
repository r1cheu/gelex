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

#include "gelex/data/simulate.h"

#include <algorithm>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader/bim_loader.h"
#include "../src/utils/formatter.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/effect_sampler.h"
#include "gelex/data/genetic_value_calculator.h"
#include "gelex/data/phenotype_generator.h"
#include "gelex/data/simulation_writer.h"
#include "gelex/exception.h"
#include "gelex/logger.h"

namespace gelex
{
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace
{
void log_simulation_params(const PhenotypeSimulator::Config& config)
{
    auto logger = logging::get();
    logger->info(section("Simulation Parameters..."));
    logger->info(
        task("Heritability (h²)      : {:.2f}", config.add_heritability));
    if (config.dom_heritability > 0.0)
    {
        logger->info(
            task("Dom-Heritability (δ²)  : {:.2f}", config.dom_heritability));
    }
    if (config.intercept != 0.0)
    {
        logger->info(task("Intercept              : {:.2f}", config.intercept));
    }
    logger->info(task("Seed                   : {}", config.seed));
    logger->info("");
}

void log_data_info(Eigen::Index n_snps, Eigen::Index n_samples)
{
    auto logger = logging::get();
    logger->info(section("Loading Data..."));
    logger->info(task("SNPs                   : {}", n_snps));
    logger->info(task("Samples                : {}", n_samples));
    logger->info("");
}

void log_phenotype_stats(double true_h2, double true_d2)
{
    auto logger = logging::get();
    logger->info(section("Generating Phenotypes..."));
    logger->info(task("True h²                : {:.4f}", true_h2));
    if (true_d2 > 0.0)
    {
        logger->info(task("True δ²                : {:.4f}", true_d2));
    }
    logger->info("");
}

}  // namespace

PhenotypeSimulator::PhenotypeSimulator(Config config)
    : config_(std::move(config))
{
    if (config_.add_heritability <= 0.0 || config_.add_heritability >= 1.0)
    {
        throw ArgumentValidationException("Heritability must be in (0, 1)");
    }
    if (config_.dom_heritability < 0.0 || config_.dom_heritability >= 1.0)
    {
        throw ArgumentValidationException(
            "Dominance variance (d2) must be in [0, 1)");
    }
    if (config_.add_heritability + config_.dom_heritability >= 1.0)
    {
        throw ArgumentValidationException("h2 + d2 must be less than 1");
    }
}

auto PhenotypeSimulator::resolve_output_path(
    const std::filesystem::path& output_path,
    const std::filesystem::path& bed_path) -> std::filesystem::path
{
    return output_path.empty() ? std::filesystem::path(bed_path) : output_path;
}

void PhenotypeSimulator::simulate()
{
    log_simulation_params(config_);

    auto bim_path = config_.bed_path;
    bim_path.replace_extension(".bim");
    detail::BimLoader bim_loader(bim_path);
    auto snp_ids = bim_loader.get_ids();

    const bool has_dominance = config_.dom_heritability > 0.0;
    EffectSampler effect_sampler({
        .add_classes = config_.add_effect_classes,
        .dom_classes = config_.dom_effect_classes,
        .has_dominance = has_dominance,
        .seed = config_.seed,
    });

    auto causal_effects
        = effect_sampler.sample(static_cast<Eigen::Index>(snp_ids.size()));

    auto sample_ptr = SampleManager::create_finalized(config_.bed_path);
    BedPipe bed_pipe(config_.bed_path, sample_ptr);

    const auto n_snps = bed_pipe.num_snps();
    const auto n_individuals = bed_pipe.num_samples();

    log_data_info(
        n_snps, static_cast<Eigen::Index>(sample_ptr->num_common_samples()));

    // Calculate genetic values using chunked processing
    VectorXd additive_values = VectorXd::Zero(n_individuals);
    VectorXd dominance_values = VectorXd::Zero(n_individuals);

    for (Eigen::Index start = 0; start < n_snps; start += SNP_CHUNK_SIZE)
    {
        const Eigen::Index end = std::min(start + SNP_CHUNK_SIZE, n_snps);

        // Check if any causal SNPs fall within this chunk
        bool has_causal_in_chunk = false;
        for (Eigen::Index i = start; i < end; ++i)
        {
            if (causal_effects.contains(i))
            {
                has_causal_in_chunk = true;
                break;
            }
        }
        if (!has_causal_in_chunk)
        {
            continue;
        }

        MatrixXd chunk = bed_pipe.load_chunk(start, end);

        auto chunk_values = GeneticValueCalculator::calculate_chunk(
            chunk, causal_effects, start, end, has_dominance);

        additive_values += chunk_values.additive;
        dominance_values += chunk_values.dominance;
    }

    GeneticValues genetic_values{
        .additive = std::move(additive_values),
        .dominance = std::move(dominance_values)};

    // Generate phenotypes
    PhenotypeGenerator phenotype_generator({
        .h2 = config_.add_heritability,
        .d2 = config_.dom_heritability,
        .intercept = config_.intercept,
        .seed = config_.seed + 1,
    });

    auto pheno_result = phenotype_generator.generate(genetic_values);

    log_phenotype_stats(pheno_result.true_h2, pheno_result.true_d2);

    // Update causal effects with the dominance scaling factor
    double dom_scale = phenotype_generator.dom_scale();
    if (dom_scale != 1.0)
    {
        for (auto& [idx, effect] : causal_effects)
        {
            effect.dominance *= dom_scale;
        }
    }

    auto output_prefix
        = resolve_output_path(config_.output_path, config_.bed_path);
    SimulationWriter writer(output_prefix);

    writer.write_phenotypes(pheno_result.phenotypes, sample_ptr->common_ids());
    writer.write_causal_effects(snp_ids, causal_effects);
}
}  // namespace gelex
