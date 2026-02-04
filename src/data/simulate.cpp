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
#include <cmath>
#include <format>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader/bim_loader.h"
#include "../src/data/parser.h"
#include "../src/utils/formatter.h"
#include "../src/utils/math_utils.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/genotype_processor.h"
#include "gelex/data/sample_manager.h"
#include "gelex/exception.h"
#include "gelex/logger.h"

namespace gelex
{
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace
{
void validate_effect_classes(
    std::span<const EffectSizeClass> classes,
    std::string_view label)
{
    if (classes.empty())
    {
        throw ArgumentValidationException(
            std::format("{} effect classes must not be empty", label));
    }

    double total_proportion = 0.0;
    for (const auto& cls : classes)
    {
        if (cls.proportion <= 0.0)
        {
            throw ArgumentValidationException(
                std::format("{} effect class proportion must be > 0", label));
        }
        if (cls.variance < 0.0)
        {
            throw ArgumentValidationException(
                std::format("{} effect class variance must be >= 0", label));
        }
        total_proportion += cls.proportion;
    }

    if (std::abs(total_proportion - 1.0) > 1e-6)
    {
        throw ArgumentValidationException(
            std::format(
                "{} effect class proportions must sum to 1.0 (got {})",
                label,
                total_proportion));
    }
}

// Partition n items into classes by proportion, assigning
// remainders to the last class, then shuffle the assignments.
auto assign_effect_classes(
    std::span<const EffectSizeClass> classes,
    Eigen::Index count,
    std::mt19937_64& rng) -> std::vector<int>
{
    std::vector<int> assignments(count);
    const auto n_classes = static_cast<int>(classes.size());
    Eigen::Index offset = 0;

    for (int cls = 0; cls < n_classes; ++cls)
    {
        const bool is_last = (cls == n_classes - 1);
        Eigen::Index class_count
            = is_last ? count - offset
                      : std::min(
                            static_cast<Eigen::Index>(std::round(
                                static_cast<double>(count)
                                * classes[cls].proportion)),
                            count - offset);

        std::fill_n(assignments.begin() + offset, class_count, cls);
        offset += class_count;
    }

    std::ranges::shuffle(assignments, rng);
    return assignments;
}

auto resolve_output_path(
    const std::filesystem::path& output_path,
    const std::filesystem::path& bed_path,
    std::string_view extension) -> std::filesystem::path
{
    auto base = output_path.empty() ? std::filesystem::path(bed_path)
                                    : std::filesystem::path(output_path);
    return base.replace_extension(extension);
}

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

void log_output_path(std::string_view label, const std::filesystem::path& path)
{
    auto logger = logging::get();
    logger->info(success(" {:<24}: {}", label, path.string()));
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

    validate_effect_classes(config_.add_effect_classes, "Additive");
    if (config_.dom_heritability > 0.0)
    {
        validate_effect_classes(config_.dom_effect_classes, "Dominance");
    }
}

void PhenotypeSimulator::simulate()
{
    initialize_rng();

    log_simulation_params(config_);

    // Read SNP IDs from .bim file
    auto bim_path = config_.bed_path;
    bim_path.replace_extension(".bim");
    detail::BimLoader bim_loader(bim_path);
    auto snp_ids = bim_loader.get_ids();

    auto causal_effects = select_causal_snps(snp_ids);

    // Setup sample manager and BedPipe
    auto fam_path = config_.bed_path;
    fam_path.replace_extension(".fam");

    SampleManager sample_manager(fam_path);
    sample_manager.finalize();
    auto sample_ptr
        = std::make_shared<SampleManager>(std::move(sample_manager));

    BedPipe bed_pipe(config_.bed_path, sample_ptr);

    log_data_info(
        static_cast<Eigen::Index>(snp_ids.size()),
        static_cast<Eigen::Index>(sample_ptr->num_common_samples()));

    auto genetic_values = calculate_genetic_values(bed_pipe, causal_effects);

    VectorXd phenotypes = generate_phenotypes(
        genetic_values.additive, genetic_values.dominance);

    write_results(phenotypes, sample_ptr);
    write_causal_effects(snp_ids, causal_effects);
}

void PhenotypeSimulator::initialize_rng()
{
    try
    {
        if (config_.seed == -1)
        {
            std::random_device rd;
            rng_.seed(rd());
        }
        else
        {
            rng_.seed(config_.seed);
        }
    }
    catch (const std::exception& e)
    {
        throw ArgumentValidationException(
            std::format("Failed to initialize RNG: {}", e.what()));
    }
}

auto PhenotypeSimulator::select_causal_snps(
    const std::vector<std::string>& snp_ids)
    -> std::unordered_map<Eigen::Index, CausalEffect>
{
    const auto n_snps = static_cast<Eigen::Index>(snp_ids.size());

    auto add_assignments
        = assign_effect_classes(config_.add_effect_classes, n_snps, rng_);

    const bool has_dominance = config_.dom_heritability > 0.0;
    auto dom_assignments
        = has_dominance
              ? assign_effect_classes(config_.dom_effect_classes, n_snps, rng_)
              : std::vector<int>{};

    // Sample effect sizes from class-specific normal distributions
    auto sample_effect
        = [&](const std::vector<EffectSizeClass>& classes, int cls) -> double
    {
        double variance = classes[cls].variance;
        if (variance == 0.0)
        {
            return 0.0;
        }
        return std::normal_distribution<double>(0.0, std::sqrt(variance))(rng_);
    };

    std::unordered_map<Eigen::Index, CausalEffect> causal_effects;
    causal_effects.reserve(n_snps);

    for (Eigen::Index i = 0; i < n_snps; ++i)
    {
        CausalEffect effect{
            .additive
            = sample_effect(config_.add_effect_classes, add_assignments[i]),
            .dominance = 0.0,
            .add_class = add_assignments[i],
            .dom_class = 0,
        };

        if (has_dominance)
        {
            effect.dominance
                = sample_effect(config_.dom_effect_classes, dom_assignments[i]);
            effect.dom_class = dom_assignments[i];
        }

        causal_effects.emplace(i, effect);
    }

    return causal_effects;
}

auto PhenotypeSimulator::calculate_genetic_values(
    BedPipe& bed_pipe,
    const std::unordered_map<Eigen::Index, CausalEffect>& causal_effects) const
    -> GeneticValues
{
    const auto n_individuals = bed_pipe.num_samples();
    const auto n_snps = bed_pipe.num_snps();
    const bool has_dominance = config_.dom_heritability > 0.0;

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

        MatrixXd dom_chunk;
        if (has_dominance)
        {
            dom_chunk = chunk;
            process_matrix<grm::OrthStandardized::Dominant>(dom_chunk);
        }

        process_matrix<grm::OrthStandardized::Additive>(chunk);

        for (const auto& [global_idx, effect] : causal_effects)
        {
            if (global_idx < start || global_idx >= end)
            {
                continue;
            }
            const Eigen::Index local_idx = global_idx - start;

            additive_values += chunk.col(local_idx) * effect.additive;

            if (has_dominance)
            {
                dominance_values += dom_chunk.col(local_idx) * effect.dominance;
            }
        }
    }

    return {
        .additive = std::move(additive_values),
        .dominance = std::move(dominance_values)};
}

auto PhenotypeSimulator::generate_phenotypes(
    const VectorXd& additive_values,
    const VectorXd& dominance_values) -> VectorXd
{
    const double h2 = config_.add_heritability;
    const double d2 = config_.dom_heritability;
    const double genetic_variance = detail::var(additive_values)(0);

    VectorXd scaled_dominance = VectorXd::Zero(additive_values.size());
    double residual_variance{};

    if (d2 > 0.0 && genetic_variance > 0.0)
    {
        // Scale dominance values so that Var(g_d) = Va * d2 / h2
        const double dominance_raw_var = detail::var(dominance_values)(0);
        const double target_dominance_var = genetic_variance * d2 / h2;

        if (dominance_raw_var > 0.0)
        {
            const double scale
                = std::sqrt(target_dominance_var / dominance_raw_var);
            scaled_dominance = dominance_values * scale;
        }

        // Ve = Va * (1 - h2 - d2) / h2
        residual_variance = genetic_variance * (1.0 - h2 - d2) / h2;
    }
    else
    {
        // d2=0: Ve = Va * (1/h2 - 1)
        residual_variance = genetic_variance * (1.0 / h2 - 1.0);
    }

    std::normal_distribution<double> residual_dist(
        0.0, std::sqrt(std::max(0.0, residual_variance)));

    VectorXd residuals(additive_values.size());
    for (Eigen::Index i = 0; i < residuals.size(); ++i)
    {
        residuals(i) = residual_dist(rng_);
    }

    VectorXd phenotypes = additive_values + scaled_dominance + residuals;
    double var_phen = detail::var(phenotypes)(0);

    true_h2_ = genetic_variance / var_phen;
    if (d2 > 0.0)
    {
        true_d2_ = detail::var(scaled_dominance)(0) / var_phen;
    }

    log_phenotype_stats(true_h2_, true_d2_);

    return phenotypes;
}

void PhenotypeSimulator::write_results(
    const VectorXd& phenotypes,
    const std::shared_ptr<SampleManager>& sample_manager) const
{
    const auto output_path
        = resolve_output_path(config_.output_path, config_.bed_path, ".phen");

    auto output = detail::open_file<std::ofstream>(output_path, std::ios::out);

    output << "FID\tIID\tphenotype\n";

    const auto& sample_ids = sample_manager->common_ids();
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(sample_ids.size());
         ++i)
    {
        const std::string_view full_id(sample_ids[i]);
        std::string_view fid = full_id;
        std::string_view iid = full_id;

        if (const auto pos = full_id.find('_'); pos != std::string_view::npos)
        {
            fid = full_id.substr(0, pos);
            iid = full_id.substr(pos + 1);
        }

        output << std::format("{}\t{}\t{}\n", fid, iid, phenotypes[i]);
    }

    log_output_path("Phenotypes saved to", output_path);
}

void PhenotypeSimulator::write_causal_effects(
    const std::vector<std::string>& snp_ids,
    const std::unordered_map<Eigen::Index, CausalEffect>& causal_effects) const
{
    auto effects_path
        = resolve_output_path(config_.output_path, config_.bed_path, ".causal");

    auto output = detail::open_file<std::ofstream>(effects_path, std::ios::out);

    output << "SNP\tadditive_effect\tdominance_effect\tadd_class\tdom_class\n";

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(snp_ids.size()); ++i)
    {
        auto it = causal_effects.find(i);
        if (it == causal_effects.end())
        {
            continue;
        }
        const auto& effect = it->second;
        output << std::format(
            "{}\t{}\t{}\t{}\t{}\n",
            snp_ids[i],
            effect.additive,
            effect.dominance,
            effect.add_class,
            effect.dom_class);
    }

    log_output_path("Causal effects saved to", effects_path);
}

}  // namespace gelex
