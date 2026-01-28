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

#include "math_utils.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace gelex
{
namespace detail
{
using Eigen::Ref;

using Eigen::Index;
using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;

RowVectorXd centralize(Ref<MatrixXd> x)
{
    RowVectorXd x_mean = RowVectorXd::Zero(x.cols());
#pragma omp parallel for default(none) shared(x, x_mean)
    for (int i = 0; i < x.cols(); ++i)
    {
        auto col = x.col(i);
        double mean_i = col.mean();
        x_mean(i) = mean_i;
        col.array() -= mean_i;
    }
    return x_mean;
}

std::pair<RowVectorXd, RowVectorXd> standardize(Ref<MatrixXd> x)
{
    RowVectorXd x_mean = RowVectorXd::Zero(x.cols());
    RowVectorXd x_stddev = RowVectorXd::Zero(x.cols());
#pragma omp parallel for default(none) shared(x, x_mean, x_stddev)
    for (int i = 0; i < x.cols(); ++i)
    {
        auto col = x.col(i);
        double mean = col.mean();
        double stddev = std::sqrt(
            (col.array() - mean).square().sum()
            / static_cast<double>(col.size() - 1));
        x_mean(i) = mean;
        x_stddev(i) = stddev;
        col.array() -= mean;
        if (stddev != 0)
        {
            col.array() /= stddev;
        }
    }
    return {x_mean, x_stddev};
}

VectorXd sum_square(const Ref<const MatrixXd>& mat)
{
    VectorXd result(mat.cols());

#pragma omp parallel for default(none) shared(result, mat)
    for (int i = 0; i < mat.cols(); ++i)
    {
        result(i) = mat.col(i).squaredNorm();
    }
    return result;
}

VectorXd sum_square(const Ref<Eigen::SparseMatrix<double>>& mat)
{
    VectorXd result(mat.cols());
#pragma omp parallel for default(none) shared(result, mat)
    for (int i = 0; i < mat.cols(); ++i)
    {
        result(i) = mat.col(i).squaredNorm();
    }
    return result;
}

}  // namespace detail
}  // namespace gelex
