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

#include "gelex/algo/sim/genetic_value_calculator.h"

#include <algorithm>
#include <cstddef>
#include <format>
#include <utility>

#include <Eigen/Core>

#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/simulate_event.h"

namespace gelex
{

namespace
{

auto notify_progress(
    const SimulateObserver& observer,
    size_t total,
    size_t current,
    bool finished = false) -> void
{
    if (!observer)
    {
        return;
    }

    SimulateEvent event;
    event.emplace<SimulateProgressEvent>(total, current, finished);
    observer(event);
}

}  // namespace

GeneticValueCalculator::GeneticValueCalculator(
    const std::filesystem::path& bed_path,
    bool has_dominance)
    : has_dominance_(has_dominance),
      sample_manager_(SampleManager::create_finalized(bed_path)),
      bed_pipe_(bed_path, sample_manager_)
{
}

auto GeneticValueCalculator::calculate(
    const CausalEffects& effects,
    const SimulateObserver& observer) const -> GeneticValues
{
    const Eigen::Index n_snps = bed_pipe_.num_snps();
    const Eigen::Index n_individuals = bed_pipe_.num_samples();

    if (effects.size() != n_snps)
    {
        throw gelex::InvalidInputException(
            std::format(
                "Number of effects should match the number of SNPs, but got {} "
                "effects for {} SNPs",
                effects.size(),
                n_snps));
    }

    Eigen::VectorXd additive_values = Eigen::VectorXd::Zero(n_individuals);
    Eigen::VectorXd dominance_values;
    if (has_dominance_)
    {
        dominance_values = Eigen::VectorXd::Zero(n_individuals);
    }
    notify_progress(observer, static_cast<size_t>(n_snps), 0);

    for (Eigen::Index start = 0; start < n_snps; start += SNP_CHUNK_SIZE)
    {
        const Eigen::Index end = std::min(start + SNP_CHUNK_SIZE, n_snps);
        auto [add_chunk, dom_chunk]
            = encode_chunk(bed_pipe_.load_chunk(start, end));
        additive_values
            += add_chunk * effects.additive.segment(start, end - start);
        if (has_dominance_)
        {
            dominance_values
                += dom_chunk * effects.dominance.segment(start, end - start);
        }
        notify_progress(
            observer, static_cast<size_t>(n_snps), static_cast<size_t>(end));
    }
    notify_progress(
        observer,
        static_cast<size_t>(n_snps),
        static_cast<size_t>(n_snps),
        true);

    return {
        .additive = std::move(additive_values),
        .dominance = std::move(dominance_values)};
}

auto GeneticValueCalculator::sample_ids() const -> std::span<const std::string>
{
    return sample_manager_->common_ids();
}

auto GeneticValueCalculator::encode_chunk(
    const Eigen::Ref<const Eigen::MatrixXd>& chunk) const
    -> std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
{
    Eigen::MatrixXd add_chunk = chunk;
    process_matrix<grm::Standardized::Additive>(add_chunk);

    Eigen::MatrixXd dom_chunk;
    if (has_dominance_)
    {
        dom_chunk = chunk;
        process_matrix<grm::Standardized::Dominant>(dom_chunk);
    }

    return {add_chunk, dom_chunk};
}

}  // namespace gelex
