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

#ifndef GELEX_INFRA_STATS_DESCRIPTIVE_H_
#define GELEX_INFRA_STATS_DESCRIPTIVE_H_

#include <type_traits>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace gelex::stats
{

[[nodiscard]] Eigen::RowVectorXd centralize(Eigen::Ref<Eigen::MatrixXd> x);
[[nodiscard]] std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd> standardize(
    Eigen::Ref<Eigen::MatrixXd> x);
[[nodiscard]] Eigen::VectorXd sum_square(
    const Eigen::Ref<const Eigen::MatrixXd>& mat);
[[nodiscard]] Eigen::VectorXd sum_square(
    const Eigen::Ref<Eigen::SparseMatrix<double>>& mat);

namespace detail
{

template <Eigen::Index Axis = 0, typename Derived>
auto var(const Eigen::DenseBase<Derived>& a, Eigen::Index norm_type = 1)
    -> std::conditional_t<Axis == 0, Eigen::RowVectorXd, Eigen::VectorXd>
{
    static_assert(Axis == 0 || Axis == 1);
    const Eigen::Index ddof = (norm_type == 0) ? 0 : 1;

    if constexpr (Axis == 0)
    {
        const Eigen::Index n = a.cols();
        Eigen::RowVectorXd result(n);
#pragma omp parallel for default(none) shared(n, a, result, ddof)
        for (Eigen::Index i = 0; i < n; ++i)
        {
            auto col = a.col(i);
            double mean_val = col.mean();
            double sum_sq = (col.array() - mean_val).square().sum();
            result(i) = sum_sq / static_cast<double>(col.size() - ddof);
        }
        return result;
    }
    else
    {
        const Eigen::Index n = a.rows();
        Eigen::VectorXd result(n);
#pragma omp parallel for default(none) shared(n, a, result, ddof)
        for (Eigen::Index i = 0; i < n; ++i)
        {
            auto row = a.row(i);
            double mean_val = row.mean();
            double sum_sq = (row.array() - mean_val).square().sum();
            result(i) = sum_sq / static_cast<double>(row.size() - ddof);
        }
        return result;
    }
}

}  // namespace detail
}  // namespace gelex::stats

#endif  // GELEX_INFRA_STATS_DESCRIPTIVE_H_
