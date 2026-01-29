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

#include <cmath>
#include <format>
#include <fstream>
#include <memory>
#include <random>
#include <sstream>
#include <string_view>
#include <vector>

#include <Eigen/Core>

#include "../src/data/parser.h"
#include "../src/utils/math_utils.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/sample_manager.h"
#include "gelex/data/variant_processor.h"
#include "gelex/exception.h"

namespace gelex
{
using Eigen::MatrixXd;
using Eigen::VectorXd;

PhenotypeSimulator::PhenotypeSimulator(Config config)
    : config_(std::move(config))
{
    if (config_.heritability <= 0.0 || config_.heritability >= 1.0)
    {
        throw ArgumentValidationException("Heritability must be in (0, 1)");
    }
}

void PhenotypeSimulator::simulate()
{
    // Stage 1: Initialization
    initialize_rng();

    auto causal_effects = load_or_generate_causal_effects();

    // Stage 2: Data Processing Pipeline Setup
    auto fam_path = config_.bed_path;
    fam_path.replace_extension(".fam");

    SampleManager sample_manager(fam_path);
    sample_manager.finalize();
    auto sample_ptr
        = std::make_shared<SampleManager>(std::move(sample_manager));

    BedPipe bed_pipe(config_.bed_path, sample_ptr);

    auto genetic_values = calculate_genetic_values(bed_pipe, causal_effects);

    VectorXd phenotypes = generate_phenotypes(genetic_values);

    write_results(phenotypes, sample_ptr);
}

// --- Private Helper Methods ---

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

auto PhenotypeSimulator::load_or_generate_causal_effects()
    -> std::unordered_map<std::string, double>
{
    auto file = detail::open_file<std::ifstream>(
        config_.causal_variants_path, std::ios::in);

    std::unordered_map<std::string, double> variants;
    std::normal_distribution<double> dist(0.0, 1.0);
    std::string line;
    bool effects_were_generated = false;

    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        std::stringstream ss(line);
        std::string id;
        double effect{};
        ss >> id;

        if (ss >> effect)
        {
            variants.emplace(std::move(id), effect);
        }
        else
        {
            effects_were_generated = true;
            variants.emplace(std::move(id), dist(rng_));
        }
    }

    // If we generated effects, write them to a new file as a side-effect.
    if (effects_were_generated)
    {
        auto effects_path = config_.causal_variants_path;
        effects_path.concat(".effects");

        auto outfile
            = detail::open_file<std::ofstream>(effects_path, std::ios::out);
        for (const auto& [id, eff] : variants)
        {
            outfile << std::format("{}\t{}\n", id, eff);
        }
    }

    return variants;
}

auto PhenotypeSimulator::calculate_genetic_values(
    BedPipe& bed_pipe,
    const std::unordered_map<std::string, double>& /* causal_effects */)
    -> VectorXd
{
    const auto n_individuals = bed_pipe.num_samples();
    const auto n_snps = bed_pipe.num_snps();

    VectorXd genetic_values = VectorXd::Zero(n_individuals);

    for (Eigen::Index start = 0; start < n_snps; start += SNP_CHUNK_SIZE)
    {
        const Eigen::Index end = std::min(start + SNP_CHUNK_SIZE, n_snps);

        MatrixXd genotype_chunk = bed_pipe.load_chunk(start, end);

        for (Eigen::Index i = 0; i < genotype_chunk.cols(); ++i)
        {
            auto col = genotype_chunk.col(i);
            const auto stats = StandardizingProcessor::process_variant(col);
            if (stats.is_monomorphic)
            {
                col.setZero();
            }
            // Note: We removed SNP ID matching since BedPipe doesn't have
            // snp_ids() This assumes all variants in the chunk are causal
            genetic_values += col;
        }
    }
    return genetic_values;
}

auto PhenotypeSimulator::generate_phenotypes(const VectorXd& genetic_values)
    -> VectorXd
{
    // Calculate variance of genetic values
    const double genetic_variance = detail::var(genetic_values)(0);

    // Derive residual variance from heritability: h^2 = Vg / (Vg + Ve)
    // => Ve = Vg * (1/h^2 - 1)
    const double residual_variance
        = genetic_variance * (1.0 / config_.heritability - 1.0);

    std::normal_distribution<double> residual_dist(
        0.0, std::sqrt(residual_variance));

    VectorXd residuals(genetic_values.size());
    for (Eigen::Index i = 0; i < residuals.size(); ++i)
    {
        residuals(i) = residual_dist(rng_);
    }

    return genetic_values + residuals;
}

void PhenotypeSimulator::write_results(
    const VectorXd& phenotypes,
    const std::shared_ptr<SampleManager>& sample_manager) const
{
    const auto output_path = config_.output_path.empty()
                                 ? std::filesystem::path(config_.bed_path)
                                       .replace_extension(".phen")
                                 : config_.output_path;

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
}

}  // namespace gelex
