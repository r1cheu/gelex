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

#ifndef GELEX_INFRA_STATS_DETAIL_VAR_H_
#define GELEX_INFRA_STATS_DETAIL_VAR_H_

#include <cassert>
#include <type_traits>

#include <Eigen/Core>

namespace gelex::stats::detail
{

template <typename Derived>
auto vecvar(const Eigen::DenseBase<Derived>& values, Eigen::Index norm_type = 1)
    -> double
{
    assert(
        (values.rows() == 1 || values.cols() == 1)
        && "vecvar: input must be a vector");

    const Eigen::Index ddof = (norm_type == 0) ? 0 : 1;
    const double mean = values.mean();
    return (values.derived().array() - mean).square().sum()
           / static_cast<double>(values.size() - ddof);
}

template <Eigen::Index Axis = 0, typename Derived>
auto matvar(const Eigen::DenseBase<Derived>& matrix, Eigen::Index norm_type = 1)
    -> std::conditional_t<Axis == 0, Eigen::RowVectorXd, Eigen::VectorXd>
{
    static_assert(Axis == 0 || Axis == 1);

    if constexpr (Axis == 0)
    {
        const Eigen::Index n = matrix.cols();
        Eigen::RowVectorXd result(n);
#pragma omp parallel for default(none) shared(n, matrix, result, norm_type)
        for (Eigen::Index i = 0; i < n; ++i)
        {
            result(i) = vecvar(matrix.col(i), norm_type);
        }
        return result;
    }
    else
    {
        const Eigen::Index n = matrix.rows();
        Eigen::VectorXd result(n);
#pragma omp parallel for default(none) shared(n, matrix, result, norm_type)
        for (Eigen::Index i = 0; i < n; ++i)
        {
            result(i) = vecvar(matrix.row(i), norm_type);
        }
        return result;
    }
}

}  // namespace gelex::stats::detail

#endif  // GELEX_INFRA_STATS_DETAIL_VAR_H_
