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

#include <cmath>
#include <cstddef>
#include <vector>

#include <nanobench.h>

#include <catch2/catch_test_macros.hpp>

#include "gelex/infra/stats/normal_sampler.h"

namespace
{

constexpr std::size_t kCount = 131'072;
using Sampler = gelex::NormalSampler<double>;

struct Inputs
{
    std::vector<double> rhs;
    std::vector<double> predictor_sum_squares;
    std::vector<double> prior_var;
    double fixed_prior_var{};
    double residual_var{};
};

auto make_inputs() -> Inputs
{
    Inputs inputs;
    inputs.rhs.resize(kCount);
    inputs.predictor_sum_squares.resize(kCount);
    inputs.prior_var.resize(kCount);
    inputs.fixed_prior_var = 0.42;
    inputs.residual_var = 1.35;

    for (std::size_t i = 0; i < kCount; ++i)
    {
        const auto x = static_cast<double>((i % 251) + 1);
        inputs.rhs[i] = (0.01 * x) - 1.2;
        inputs.predictor_sum_squares[i] = 30.0 + static_cast<double>(i % 997);
        inputs.prior_var[i] = 0.05 + (0.001 * static_cast<double>(i % 389));
    }

    return inputs;
}

}  // namespace

TEST_CASE("NormalSampler posterior throughput", "[!benchmark][stats][normal]")
{
    const auto inputs = make_inputs();

    ankerl::nanobench::Bench b;
    b.title("NormalSampler posterior")
        .unit("marker")
        .batch(kCount)
        .warmup(5)
        .minEpochIterations(20);

    b.run(
        "fixed_prior_var",
        [&]()
        {
            double sink = 0.0;
            Sampler sampler{inputs.fixed_prior_var};
            for (std::size_t i = 0; i < kCount; ++i)
            {
                const auto posterior = sampler.posterior_with_logL({
                    .quadratic = inputs.predictor_sum_squares[i],
                    .linear = inputs.rhs[i],
                    .scale = inputs.residual_var,
                });
                sink += posterior.params.mean + std::sqrt(posterior.params.var)
                        + posterior.log_likelihood_kernel;
            }
            ankerl::nanobench::doNotOptimizeAway(sink);
        });

    b.run(
        "varying_prior_var",
        [&]()
        {
            double sink = 0.0;
            Sampler sampler{inputs.prior_var[0]};
            for (std::size_t i = 0; i < kCount; ++i)
            {
                sampler.set_prior_var(inputs.prior_var[i]);
                const auto posterior = sampler.posterior_with_logL({
                    .quadratic = inputs.predictor_sum_squares[i],
                    .linear = inputs.rhs[i],
                    .scale = inputs.residual_var,
                });
                sink += posterior.params.mean + std::sqrt(posterior.params.var)
                        + posterior.log_likelihood_kernel;
            }
            ankerl::nanobench::doNotOptimizeAway(sink);
        });
}
