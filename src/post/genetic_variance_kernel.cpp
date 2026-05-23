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

#include "gelex/post/genetic_variance_kernel.h"

#include <cstddef>
#include <ranges>
#include <string>
#include <string_view>

#include <fmt/format.h>

#include "gelex/exception.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/post/detail/utils.h"

namespace gelex
{

namespace
{

struct GeneticBlockSlot
{
    std::size_t block_index;
    std::size_t slot;
};

auto mode_name(GeneticMode mode) -> std::string_view
{
    for (const auto& [value, name] : kGeneticModeNames)
    {
        if (value == mode)
        {
            return name;
        }
    }
    return "unknown";
}

auto find_genetic_index(
    const io::detail::BinaryReader& reader,
    GeneticMode kind) -> std::size_t
{
    constexpr std::string_view modes_path = "genetic/modes";
    if (!reader.contains(modes_path))
    {
        throw GelexException("samples do not contain genetic/modes");
    }

    const auto target = mode_name(kind);
    const auto modes = reader.to_strings(modes_path);
    for (const auto [index, mode] : std::views::enumerate(modes))
    {
        if (mode == target)
        {
            return static_cast<std::size_t>(index);
        }
    }

    throw GelexException(
        fmt::format("samples do not contain genetic mode {}", kind));
}

auto find_genetic_block_slot(
    const io::detail::BinaryReader& reader,
    GeneticMode kind) -> GeneticBlockSlot
{
    const auto target = mode_name(kind);
    for (std::size_t block_index = 0;; ++block_index)
    {
        const auto modes_path
            = fmt::format("genetic_block/{}/modes", block_index);
        if (!reader.contains(modes_path))
        {
            break;
        }

        const auto modes = reader.to_strings(modes_path);
        for (const auto [slot, mode] : std::views::enumerate(modes))
        {
            if (mode == target)
            {
                return {
                    .block_index = block_index,
                    .slot = static_cast<std::size_t>(slot)};
            }
        }
    }

    throw GelexException(
        fmt::format("samples do not contain genetic block for {}", kind));
}

}  // namespace

GeneticVarianceProcessor::GeneticVarianceProcessor(
    const GeneticInput& input,
    std::span<const io::detail::BinaryReader> readers)
    : matrix_{input.genotype->matrix()}, kind_{input.kind}, n_components_{0}
{
    const auto& ref = readers.front();
    const auto genetic_index = find_genetic_index(ref, kind_);
    const auto block_slot = find_genetic_block_slot(ref, kind_);
    const auto coeff_path = fmt::format("genetic/{}/coeffs", genetic_index);
    const auto tracker_path = fmt::format(
        "genetic_block/{}/prior_state/proportion/{}/assignment",
        block_slot.block_index,
        block_slot.slot);

    if (ref.contains(tracker_path))
    {
        n_components_ = static_cast<Eigen::Index>(
            ref.to_map<int>(tracker_path).maxCoeff());
    }

    auto n_records = ref.to_map<double>(coeff_path).cols();

    for (const auto& reader : readers)
    {
        coeff_maps_.push_back(reader.to_map<double>(coeff_path));
        if (n_components_ > 0)
        {
            tracker_maps_.push_back(reader.to_map<int>(tracker_path));
            component_variance_chains_.emplace_back(n_components_, n_records);
        }
    }
}

auto GeneticVarianceProcessor::process(
    std::size_t chain_idx,
    Eigen::Index col_begin,
    Eigen::Index col_end) -> const Eigen::MatrixXd&
{
    const auto chunk_cols = col_end - col_begin;
    const auto beta_block
        = coeff_maps_[chain_idx].middleCols(col_begin, chunk_cols);

    if (n_components_ == 0)
    {
        gebv_chunk_.noalias() = matrix_ * beta_block;
        last_variances_ = stats::detail::var(gebv_chunk_);
        return gebv_chunk_;
    }

    gebv_chunk_.setZero(matrix_.rows(), chunk_cols);
    const auto tracker_block
        = tracker_maps_[chain_idx].middleCols(col_begin, chunk_cols);
    masked_beta_chunk_.resize(beta_block.rows(), chunk_cols);

    for (Eigen::Index k = 1; k <= n_components_; ++k)
    {
        // select beta where tracker == k, else 0
        masked_beta_chunk_ = (tracker_block.array() == k)
                                 .cast<double>()
                                 .cwiseProduct(beta_block.array())
                                 .matrix();
        component_gebv_chunk_.noalias() = matrix_ * masked_beta_chunk_;
        component_variance_chains_[chain_idx]
            .middleCols(col_begin, chunk_cols)
            .row(k - 1)
            = stats::detail::var(component_gebv_chunk_);
        gebv_chunk_ += component_gebv_chunk_;
    }

    last_variances_ = stats::detail::var(gebv_chunk_);
    return gebv_chunk_;
}

auto GeneticVarianceProcessor::build_diagnostics(double hdpi_threshold) const
    -> std::vector<ParameterDiag>
{
    if (n_components_ == 0)
    {
        return {};
    }

    auto diags = post::detail::compute_posterior_summaries(
        component_variance_chains_, hdpi_threshold);

    for (Eigen::Index k = 0; k < n_components_; ++k)
    {
        diags[static_cast<size_t>(k)].section = fmt::format("{}", kind_);
        diags[static_cast<size_t>(k)].name = fmt::format("σ²[{}]", k + 1);
    }
    return diags;
}

}  // namespace gelex
