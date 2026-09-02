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
#include "gelex/data/rank_inverse_norm_transform.h"

#include <Eigen/Core>
#include <algorithm>
#include <numeric>
#include <ranges>
#include <vector>

#include "gelex/infra/normal.h"

namespace gelex
{

namespace
{

auto compute_ranks(const Eigen::Ref<const Eigen::VectorXd>& values)
    -> Eigen::VectorXd
{
    Eigen::Index n = values.size();

    std::vector<Eigen::Index> indices(n);
    std::ranges::iota(indices, Eigen::Index{0});
    std::ranges::sort(
        indices,
        [&](Eigen::Index a, Eigen::Index b) { return values[a] < values[b]; });

    Eigen::VectorXd ranks(n);

    Eigen::Index i = 0;
    while (i < n)
    {
        Eigen::Index j = i;
        while (j < n && values[indices[j]] == values[indices[i]])
        {
            ++j;
        }

        double avg_rank = static_cast<double>(i + 1 + j) / 2.0;
        for (Eigen::Index k = i; k < j; ++k)
        {
            ranks[indices[k]] = avg_rank;
        }

        i = j;
    }

    return ranks;
}

auto compute_residuals(
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::MatrixXd>& X) -> Eigen::VectorXd
{
    return y - X * (X.transpose() * X).ldlt().solve(X.transpose() * y);
}

}  // namespace

auto direct_int(Eigen::Ref<Eigen::VectorXd> phenotype, double offset) -> void
{
    auto ranks = compute_ranks(phenotype);
    auto n = static_cast<double>(phenotype.size());
    double denominator = n - (2.0 * offset) + 1.0;

    for (auto&& [rank, value] : std::views::zip(ranks, phenotype))
    {
        double quantile
            = std::clamp((rank - offset) / denominator, 1e-10, 1.0 - 1e-10);
        value = norm_ppf(quantile);
    }
}

auto indirect_int(
    Eigen::Ref<Eigen::VectorXd> phenotype,
    const Eigen::Ref<const Eigen::MatrixXd>& covariates,
    double offset) -> void
{
    auto residuals = compute_residuals(phenotype, covariates);
    direct_int(residuals, offset);
    phenotype = residuals;
}

}  // namespace gelex
