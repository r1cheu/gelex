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

#include <nanobench.h>

#include <Eigen/Dense>

int main()
{
    const Eigen::Index M = 2000;
    const Eigen::Index N = 10000;

    Eigen::MatrixXd Z = Eigen::MatrixXd::Random(M, N);
    Eigen::MatrixXd V_inv = Eigen::MatrixXd::Random(M, M);
    V_inv = V_inv * V_inv.transpose();
    Eigen::MatrixXd V_inv_Z = V_inv * Z;

    ankerl::nanobench::Bench b;

    // 方案 0: Naive Full Multiply
    b.run(
        "Zt_Vinv_Z M=2000 N=10000 naive_full_multiply",
        [&]()
        {
            Eigen::VectorXd result = (Z.transpose() * V_inv_Z).diagonal();
            ankerl::nanobench::doNotOptimizeAway(result);
        });

    // 方案 1: Row-wise Dot (Eigen Optimized)
    b.run(
        "Zt_Vinv_Z M=2000 N=10000 row_wise_dot",
        [&]()
        {
            Eigen::VectorXd result(N);
            for (Eigen::Index i = 0; i < N; ++i)
            {
                result(i) = V_inv_Z.row(i).dot(Z.col(i));
            }
            ankerl::nanobench::doNotOptimizeAway(result);
        });

    // 方案 2: Col-wise Loop
    b.run(
        "Zt_Vinv_Z M=2000 N=10000 col_wise_loop",
        [&]()
        {
            Eigen::VectorXd result(N);
            for (Eigen::Index i = 0; i < N; ++i)
            {
                result(i) = Z.col(i).dot(V_inv * Z.col(i));
            }
            ankerl::nanobench::doNotOptimizeAway(result);
        });

    // 方案 3: Array colwise sum
    b.run(
        "Zt_Vinv_Z M=2000 N=10000 array_colwise_sum",
        [&]()
        {
            Eigen::VectorXd result
                = (Z.array() * V_inv_Z.array()).colwise().sum();
            ankerl::nanobench::doNotOptimizeAway(result);
        });
}
