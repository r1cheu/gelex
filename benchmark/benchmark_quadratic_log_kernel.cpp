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

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <nanobench.h>
#include <vector>

#include "gelex/bayes/stats/quadratic_log_kernel.h"

namespace
{

constexpr std::size_t marker_count = 131'072;

struct Inputs
{
    std::vector<double> rhs;
    std::vector<double> predictor_sum_squares;
    std::vector<double> prior_variances;
    double fixed_prior_variance{};
    double residual_variance{};
};

auto make_inputs() -> Inputs
{
    Inputs inputs;
    inputs.rhs.resize(marker_count);
    inputs.predictor_sum_squares.resize(marker_count);
    inputs.prior_variances.resize(marker_count);
    inputs.fixed_prior_variance = 0.42;
    inputs.residual_variance = 1.35;

    for (std::size_t i = 0; i < marker_count; ++i)
    {
        const auto x = static_cast<double>((i % 251) + 1);
        inputs.rhs[i] = (0.01 * x) - 1.2;
        inputs.predictor_sum_squares[i] = 30.0 + static_cast<double>(i % 997);
        inputs.prior_variances[i]
            = 0.05 + (0.001 * static_cast<double>(i % 389));
    }

    return inputs;
}

}  // namespace

TEST_CASE(
    "QuadraticLogKernel posterior throughput",
    "[!benchmark][stats][quadratic_log_kernel]")
{
    const auto inputs = make_inputs();

    ankerl::nanobench::Bench b;
    b.title("QuadraticLogKernel posterior")
        .unit("marker")
        .batch(marker_count)
        .warmup(5)
        .minEpochIterations(20);

    b.run(
        "fixed_prior_variance",
        [&]()
        {
            double sink = 0.0;
            const auto prior
                = gelex::make_normal_prior(inputs.fixed_prior_variance);
            for (std::size_t i = 0; i < marker_count; ++i)
            {
                const auto likelihood = gelex::make_coefficient_likelihood(
                    inputs.predictor_sum_squares[i],
                    inputs.rhs[i],
                    inputs.residual_variance);
                const auto posterior = likelihood + prior;
                const auto parameters = posterior.normal_parameters();
                sink += parameters.mean() + parameters.stddev()
                        + posterior.log_integral();
            }
            ankerl::nanobench::doNotOptimizeAway(sink);
        });

    b.run(
        "varying_prior_variance",
        [&]()
        {
            double sink = 0.0;
            for (std::size_t i = 0; i < marker_count; ++i)
            {
                const auto likelihood = gelex::make_coefficient_likelihood(
                    inputs.predictor_sum_squares[i],
                    inputs.rhs[i],
                    inputs.residual_variance);
                const auto posterior
                    = likelihood
                      + gelex::make_normal_prior(inputs.prior_variances[i]);
                const auto parameters = posterior.normal_parameters();
                sink += parameters.mean() + parameters.stddev()
                        + posterior.log_integral();
            }
            ankerl::nanobench::doNotOptimizeAway(sink);
        });
}
