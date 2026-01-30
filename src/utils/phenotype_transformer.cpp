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
#include "phenotype_transformer.h"

#include <algorithm>
#include <numeric>
#include <ranges>
#include <vector>

#include "math_utils.h"

namespace gelex::detail
{

PhenotypeTransformer::PhenotypeTransformer(double offset) : offset_(offset) {}

auto PhenotypeTransformer::apply_dint(
    Eigen::Ref<Eigen::VectorXd> phenotype) const -> void
{
    int_transform(phenotype);
}

auto PhenotypeTransformer::apply_iint(
    Eigen::Ref<Eigen::VectorXd> phenotype,
    const Eigen::Ref<const Eigen::MatrixXd>& covariates) const -> void
{
    auto residuals = compute_residuals(phenotype, covariates);

    int_transform(residuals);
    phenotype = residuals;
}

auto PhenotypeTransformer::int_transform(
    Eigen::Ref<Eigen::VectorXd> values) const -> void
{
    auto ranks = compute_ranks(values);
    double n = static_cast<double>(values.size());
    double denominator = n - 2.0 * offset_ + 1.0;

    for (auto&& [rank, value] : std::views::zip(ranks, values))
    {
        double quantile
            = std::clamp((rank - offset_) / denominator, 1e-10, 1.0 - 1e-10);
        value = inverse_of_normal_cdf(quantile, 0.0, 1.0);
    }
}

auto PhenotypeTransformer::compute_ranks(
    const Eigen::Ref<const Eigen::VectorXd>& values) -> Eigen::VectorXd
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

        double avg_rank = (i + 1 + j) / 2.0;
        for (Eigen::Index k = i; k < j; ++k)
        {
            ranks[indices[k]] = avg_rank;
        }

        i = j;
    }

    return ranks;
}

auto PhenotypeTransformer::compute_residuals(
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::MatrixXd>& X) -> Eigen::VectorXd
{
    return y - X * (X.transpose() * X).ldlt().solve(X.transpose() * y);
}

}  // namespace gelex::detail
