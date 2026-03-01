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

#include "gelex/pipeline/post/posterior_analysis_engine.h"

#include <fmt/format.h>

#include <Eigen/Core>
#include <cmath>
#include <string>
#include <vector>

#include "gelex/algo/stats/diagnostics.h"
#include "gelex/data/io/binary_mmap_loader.h"
#include "gelex/exception.h"  // FileFormatException, InvalidInputException
#include "gelex/infra/logging/post_event.h"

namespace
{

constexpr Eigen::Index kIdxResidVar = 0;
constexpr Eigen::Index kIdxAddVar = 1;
constexpr Eigen::Index kIdxAddH2 = 2;
constexpr Eigen::Index kIdxDomVar = 3;
constexpr Eigen::Index kIdxDomH2 = 4;

auto compute_diag(
    Eigen::Index row,
    std::string_view name,
    const gelex::Chains& chains,
    const Eigen::VectorXd& ess_vec,
    const Eigen::MatrixXd& rhat_mat) -> gelex::ParamDiag
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
    auto [lo, hi] = gelex::hpdi(all, 0.9);

    return gelex::ParamDiag{
        .name = std::string(name),
        .mean = mean,
        .sd = sd,
        .hpdi_lo = lo,
        .hpdi_hi = hi,
        .ess = ess_vec(row),
        .rhat = rhat_mat(row, 0)};
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

    for (const auto& prefix : config_.in_prefixes)
    {
        auto path = prefix + ".scalar_chain";
        detail::BinaryMmapLoader<double> loader(path);
        const auto& mat = loader.matrix();

        if (n_scalars < 0)
        {
            n_scalars = mat.rows();
            n_records = mat.cols();

            if (n_scalars != 3 && n_scalars != 5)
            {
                throw FileFormatException(
                    fmt::format(
                        "Unexpected scalar_chain shape ({} rows) in '{}'",
                        n_scalars,
                        path));
            }
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

        chains.emplace_back(loader.load_copy());
    }

    const bool has_dominance = (n_scalars == 5);

    auto ess = effect_sample_size(chains);
    auto rhat = split_gelman_rubin(chains);

    std::vector<ParamDiag> diags;
    diags.reserve(static_cast<size_t>(n_scalars));

    diags.push_back(compute_diag(kIdxResidVar, "σ²_e", chains, ess, rhat));
    diags.push_back(compute_diag(kIdxAddVar, "σ²_add", chains, ess, rhat));
    diags.push_back(compute_diag(kIdxAddH2, "h²_add", chains, ess, rhat));

    if (has_dominance)
    {
        diags.push_back(compute_diag(kIdxDomVar, "σ²_dom", chains, ess, rhat));
        diags.push_back(compute_diag(kIdxDomH2, "h²_dom", chains, ess, rhat));
    }

    if (observer)
    {
        observer(
            DiagnosticsReadyEvent{
                .diags = std::move(diags),
                .n_chains = n_chains,
                .n_records = n_records});
    }
}

}  // namespace gelex
