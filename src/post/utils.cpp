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

#include "gelex/post/utils.h"

#include <Eigen/Core>

#include "gelex/infra/logging/post_event.h"
#include "gelex/infra/stats/diagnostics.h"

namespace gelex::post::detail
{

auto assemble_chains(
    std::span<const ::gelex::io::detail::BinaryReader> readers,
    std::string_view section_path) -> stats::Chains
{
    stats::Chains chains;
    chains.reserve(readers.size());
    for (const auto& reader : readers)
    {
        chains.emplace_back(reader.to_mat<double>(section_path));
    }
    return chains;
}

auto compute_posterior_summaries(
    const stats::Chains& samples,
    double hdpi_threshold) -> std::vector<ParameterDiag>
{
    auto ess = stats::effect_sample_size(samples);
    auto r_hat = stats::split_gelman_rubin(samples);
    auto [intervals, medians] = stats::hpdi(samples, hdpi_threshold);

    Eigen::Index n_params = samples[0].rows();
    Eigen::Index total = 0;
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(n_params);
    Eigen::VectorXd sum_sq = Eigen::VectorXd::Zero(n_params);

    for (const auto& chain : samples)
    {
        sum += chain.rowwise().sum();
        sum_sq += chain.array().square().matrix().rowwise().sum();
        total += chain.cols();
    }

    Eigen::VectorXd mean = sum / static_cast<double>(total);
    Eigen::VectorXd sd = ((sum_sq.array() - mean.array().square() * total)
                          / static_cast<double>(total - 1))
                             .sqrt();

    std::vector<ParameterDiag> diags(static_cast<size_t>(n_params));
    for (Eigen::Index i{0}; i < n_params; ++i)
    {
        diags[static_cast<size_t>(i)].mean = mean(i);
        diags[static_cast<size_t>(i)].median = medians(i);
        diags[static_cast<size_t>(i)].sd = sd(i);
        diags[static_cast<size_t>(i)].hpdi_lo = intervals(i, 0);
        diags[static_cast<size_t>(i)].hpdi_hi = intervals(i, 1);
        diags[static_cast<size_t>(i)].ess = ess(i);
        diags[static_cast<size_t>(i)].rhat = r_hat(i);
    }
    return diags;
}

auto summarize_section(
    std::span<const ::gelex::io::detail::BinaryReader> readers,
    std::string_view section_path,
    double hdpi_threshold,
    std::string_view section,
    std::string_view name) -> ParameterDiag
{
    auto chains = assemble_chains(readers, section_path);
    auto diags = compute_posterior_summaries(chains, hdpi_threshold);
    diags[0].section = section;
    diags[0].name = name;
    return std::move(diags[0]);
}

auto summarize_section(
    std::span<const ::gelex::io::detail::BinaryReader> readers,
    std::string_view section_path,
    double hdpi_threshold,
    std::string_view section,
    std::span<const std::string> names) -> std::vector<ParameterDiag>
{
    auto chains = assemble_chains(readers, section_path);
    auto diags = compute_posterior_summaries(chains, hdpi_threshold);
    for (size_t i{0}; i < names.size(); ++i)
    {
        diags[i].section = section;
        diags[i].name = names[i];
    }
    return diags;
}

}  // namespace gelex::post::detail
