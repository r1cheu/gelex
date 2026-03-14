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

#include "gelex/pipeline/posterior_analysis_engine.h"

#include <cmath>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <Eigen/Core>

#include "gelex/algo/stats/diagnostics.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/logging/post_event.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_reader.h"

namespace
{

using gelex::Chains;
using gelex::EffectType;
using gelex::hpdi;
using gelex::ParamDiag;
using gelex::detail::BinaryReader;
using gelex::detail::DataKind;

auto compute_diag(
    Eigen::Index row,
    std::string_view name,
    const Chains& chains,
    const Eigen::VectorXd& ess_vec,
    const Eigen::MatrixXd& rhat_mat) -> ParamDiag
{
    const auto n_chains = static_cast<Eigen::Index>(chains.size());
    const Eigen::Index n_records = chains[0].cols();

    Eigen::VectorXd all(n_chains * n_records);
    for (Eigen::Index c = 0; c < n_chains; ++c)
    {
        all.segment(c * n_records, n_records) = chains[c].row(row).transpose();
    }

    double mean = all.mean();
    double sd = std::sqrt((all.array() - mean).square().mean());
    auto [lo, hi] = hpdi(all, 0.9);

    return ParamDiag{
        .name = std::string(name),
        .mean = mean,
        .sd = sd,
        .hpdi_lo = lo,
        .hpdi_hi = hi,
        .ess = ess_vec(row),
        .rhat = rhat_mat(row, 0)};
}

struct ScalarChainResult
{
    Eigen::MatrixXd data;
    Eigen::Index n_random;
    bool has_add;
    bool has_dom;
};

// Build scalar chains from BinaryReader sections.
// Row layout: residual_var, [random_0_var, ...], [add_var, add_h2], [dom_var,
// dom_h2] h2 = genetic_var / total_var  (computed on-the-fly)
auto build_scalar_chain(const BinaryReader& reader) -> ScalarChainResult
{
    // Residual variance (zero-copy mmap view)
    auto resid_var
        = reader.map<double>(EffectType::residual(), DataKind::Variance);
    const Eigen::Index n_records = resid_var.cols();

    // Count random effect variances
    Eigen::Index n_random = 0;
    for (uint8_t i = 0;; ++i)
    {
        if (!reader.has_section(EffectType::random(), DataKind::Variance, i))
        {
            break;
        }
        ++n_random;
    }

    const bool has_add
        = reader.has_section(EffectType::add(), DataKind::Variance);
    const bool has_dom
        = reader.has_section(EffectType::dom(), DataKind::Variance);

    // Count rows: 1 (resid) + n_random + 2 per genetic (var + h2)
    Eigen::Index n_rows = 1 + n_random + (has_add ? 2 : 0) + (has_dom ? 2 : 0);

    Eigen::MatrixXd result(n_rows, n_records);

    // Compute total variance per sample for h2 calculation
    Eigen::RowVectorXd total_var = resid_var.row(0);
    for (Eigen::Index i = 0; i < n_random; ++i)
    {
        total_var += reader
                         .map<double>(
                             EffectType::random(),
                             DataKind::Variance,
                             static_cast<uint8_t>(i))
                         .row(0);
    }
    if (has_add)
    {
        total_var
            += reader.map<double>(EffectType::add(), DataKind::Variance).row(0);
    }
    if (has_dom)
    {
        total_var
            += reader.map<double>(EffectType::dom(), DataKind::Variance).row(0);
    }

    // Fill result matrix
    Eigen::Index row = 0;
    result.row(row++) = resid_var.row(0);

    for (Eigen::Index i = 0; i < n_random; ++i)
    {
        result.row(row++) = reader
                                .map<double>(
                                    EffectType::random(),
                                    DataKind::Variance,
                                    static_cast<uint8_t>(i))
                                .row(0);
    }

    if (has_add)
    {
        auto add_var
            = reader.map<double>(EffectType::add(), DataKind::Variance);
        result.row(row++) = add_var.row(0);
        result.row(row++) = add_var.row(0).array() / total_var.array();
    }
    if (has_dom)
    {
        auto dom_var
            = reader.map<double>(EffectType::dom(), DataKind::Variance);
        result.row(row++) = dom_var.row(0);
        result.row(row++) = dom_var.row(0).array() / total_var.array();
    }

    return {std::move(result), n_random, has_add, has_dom};
}

}  // namespace

namespace gelex
{

PosteriorAnalysisEngine::PosteriorAnalysisEngine(Config config)
    : config_(std::move(config))
{
}

auto PosteriorAnalysisEngine::run(const PostObserver& observer) -> void
{
    const auto n_chains = static_cast<Eigen::Index>(config_.in_prefixes.size());

    Chains chains;
    chains.reserve(static_cast<size_t>(n_chains));

    Eigen::Index n_scalars = -1;
    Eigen::Index n_records = -1;

    Eigen::Index n_random = 0;
    bool has_add = false;
    bool has_dom = false;

    for (const auto& prefix : config_.in_prefixes)
    {
        auto path = prefix + ".samples";
        detail::BinaryReader reader(path);
        auto [mat, nr, ha, hd] = build_scalar_chain(reader);

        if (n_scalars < 0)
        {
            n_scalars = mat.rows();
            n_records = mat.cols();
            n_random = nr;
            has_add = ha;
            has_dom = hd;
        }
        else if (mat.rows() != n_scalars || mat.cols() != n_records)
        {
            throw InvalidInputException(
                fmt::format(
                    "Scalar chain shape mismatch in '{}': expected {}x{}, got "
                    "{}x{}",
                    path,
                    n_scalars,
                    n_records,
                    mat.rows(),
                    mat.cols()));
        }

        chains.emplace_back(std::move(mat));
    }

    auto ess = effect_sample_size(chains);
    auto rhat = split_gelman_rubin(chains);

    std::vector<ParamDiag> diags;
    diags.reserve(static_cast<size_t>(n_scalars));

    Eigen::Index row = 0;
    diags.push_back(compute_diag(row++, "σ²_e", chains, ess, rhat));

    // Skip random variance rows (no named diagnostics for them)
    row += n_random;

    if (has_add)
    {
        diags.push_back(compute_diag(row++, "σ²_add", chains, ess, rhat));
        diags.push_back(compute_diag(row++, "h²_add", chains, ess, rhat));
    }

    if (has_dom)
    {
        diags.push_back(compute_diag(row++, "σ²_dom", chains, ess, rhat));
        diags.push_back(compute_diag(row++, "h²_dom", chains, ess, rhat));
    }

    notify(
        observer,
        DiagnosticsReadyEvent{
            .diags = std::move(diags),
            .n_chains = n_chains,
            .n_records = n_records});
}

}  // namespace gelex
