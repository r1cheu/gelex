// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_INFRA_VAR_H_
#define GELEX_INFRA_VAR_H_

#include <Eigen/Core>
#include <cassert>
#include <cstdint>
#include <type_traits>

namespace gelex
{
enum class VarNormType : std::uint8_t
{
    Population,
    Sample
};

template <typename Derived>
auto vecvar(
    const Eigen::DenseBase<Derived>& values,
    VarNormType norm_type = VarNormType::Sample) -> double
{
    assert(
        (values.rows() == 1 || values.cols() == 1)
        && "vecvar: input must be a vector");

    const Eigen::Index ddof = (norm_type == VarNormType::Population) ? 0 : 1;
    const double mean = values.mean();
    return (values.derived().array() - mean).square().sum()
           / static_cast<double>(values.size() - ddof);
}

template <Eigen::Index Axis = 0, typename Derived>
auto matvar(
    const Eigen::DenseBase<Derived>& matrix,
    VarNormType norm_type = VarNormType::Sample)
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

}  // namespace gelex

#endif  // GELEX_INFRA_VAR_H_
