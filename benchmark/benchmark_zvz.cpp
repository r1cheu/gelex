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

#include <benchmark/benchmark.h>
#include <Eigen/Dense>

static void BM_Zt_Vinv_Z_Diagonal(benchmark::State& state)
{
    const Eigen::Index M = state.range(0);
    const Eigen::Index N = state.range(1);

    Eigen::MatrixXd Z = Eigen::MatrixXd::Random(M, N);
    Eigen::MatrixXd V_inv = Eigen::MatrixXd::Random(M, M);
    V_inv = V_inv * V_inv.transpose();  // 确保对称正定
    //
    Eigen::MatrixXd V_inv_Z = V_inv * Z;

    for (auto _ : state)
    {
        Eigen::VectorXd result;
        if (state.range(2) == 0)
        {  // 方案 0: Naive Full Multiply
            benchmark::DoNotOptimize(
                result = (Z.transpose() * V_inv_Z).diagonal());
        }
        else if (state.range(2) == 1)
        {  // 方案 1: Row-wise Dot (Eigen Optimized)
            // 先算 Z^T * V_inv -> (N x M)
            // 然后与 Z^T -> (N x M) 进行行对行的点积
            result.resize(N);
            for (Eigen::Index i = 0; i < N; ++i)
            {
                result(i) = V_inv_Z.row(i).dot(Z.col(i));
            }
            benchmark::DoNotOptimize(result);
        }
        else if (state.range(2) == 2)
        {  // 方案 2: Col-wise Loop
            //
            result.resize(N);
            for (Eigen::Index i = 0; i < N; ++i)
            {
                // z_i^T * (V_inv * z_i)
                result(i) = Z.col(i).dot(V_inv * Z.col(i));
            }
            benchmark::DoNotOptimize(result);
        }
        else if (state.range(2) == 3)
        {  // 方案 2: Col-wise Loop
            //
            result.resize(N);
            benchmark::DoNotOptimize(
                result = (Z.array() * V_inv_Z.array()).colwise().sum());
        }
    }
}

// 注册测试：M=2000 样本, N=500 SNPs
// 最后一个参数：0=Naive, 1=Row-wise, 2=Col-wise
BENCHMARK(BM_Zt_Vinv_Z_Diagonal)
    ->Args({2000, 10000, 0})
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Zt_Vinv_Z_Diagonal)
    ->Args({2000, 10000, 1})
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Zt_Vinv_Z_Diagonal)
    ->Args({2000, 10000, 2})
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Zt_Vinv_Z_Diagonal)
    ->Args({2000, 10000, 3})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
