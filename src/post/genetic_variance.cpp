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

#include "gelex/post/genetic_variance.h"

#include <ranges>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/algorithm_shape.h"
#include "gelex/post/detail/utils.h"

namespace gelex
{

GeneticVariancePosteriorProcessor::GeneticVariancePosteriorProcessor(
    std::span<const io::detail::BinaryReader> readers,
    std::span<const GeneticInput> genetics,
    double hdpi_threshold)
    : readers_{readers}, genetics_{genetics}, hdpi_threshold_{hdpi_threshold}
{
}

auto GeneticVariancePosteriorProcessor::process() -> GebvVarianceResult
{
    std::vector<GeneticVarianceProcessor> processors;
    processors.reserve(genetics_.size());
    for (const auto& input : genetics_)
    {
        processors.emplace_back(input, readers_);
    }

    const auto n_individuals = genetics_.front().genotype->rows();
    const auto n_chains = readers_.size();
    const auto& ref = readers_.front();
    const auto n_records
        = ref.to_map<double>(fmt::format("{}/coeff", genetics_.front().kind))
              .cols();

    constexpr Eigen::Index kChunkSize = 64;

    const auto n_variances = static_cast<Eigen::Index>(genetics_.size()) + 1;
    stats::Chains genetic_variances(
        n_chains, Eigen::MatrixXd(n_variances, n_records));

    for (std::size_t ci = 0; ci < n_chains; ++ci)
    {
        auto& variance = genetic_variances[ci];

        for (Eigen::Index col_begin = 0; col_begin < n_records;
             col_begin += kChunkSize)
        {
            const auto col_end = std::min(col_begin + kChunkSize, n_records);
            const auto chunk_cols = col_end - col_begin;

            Eigen::MatrixXd gebv_total
                = Eigen::MatrixXd::Zero(n_individuals, chunk_cols);

            for (auto&& [kind_idx, processor] :
                 std::views::enumerate(processors))
            {
                gebv_total += processor.process(ci, col_begin, col_end);
                variance.row(static_cast<Eigen::Index>(kind_idx))
                    .segment(col_begin, chunk_cols)
                    = processor.last_variances().head(chunk_cols);
            }

            variance.row(n_variances - 1).segment(col_begin, chunk_cols)
                = stats::detail::var(gebv_total).head(chunk_cols);
        }
    }

    // diagnostics
    std::vector<ParameterDiag> diags;
    for (const auto& processor : processors)
    {
        diags.append_range(processor.build_diagnostics(hdpi_threshold_));
    }

    auto total_diags = post::detail::compute_posterior_summaries(
        genetic_variances, hdpi_threshold_);

    for (auto&& [kind_idx, input] : std::views::enumerate(genetics_))
    {
        auto& diag = total_diags[static_cast<size_t>(kind_idx)];
        diag.section = "Genetic";
        diag.name = bayes::to_variance_label(input.kind);
    }
    total_diags.back().section = "Genetic";
    total_diags.back().name = "σ²_total";

    diags.append_range(std::move(total_diags));

    return {std::move(diags), std::move(genetic_variances)};
}

}  // namespace gelex
