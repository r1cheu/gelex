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

#include <Eigen/Core>

#include "gelex/algo/infer/detail/marker_op.h"
#include "gelex/infra/stats/conjugate_prior.h"

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

void require_equivalent_outputs(const Inputs& inputs)
{
    Sampler fixed_sampler{inputs.fixed_prior_var};
    Sampler varying_sampler{inputs.prior_var[0]};
    const double fixed_residual_over_prior_var
        = inputs.residual_var / inputs.fixed_prior_var;
    bool outputs_equivalent = true;

    for (std::size_t i = 0; i < 1024; ++i)
    {
        const Sampler::Kernel kernel{
            .quadratic = inputs.predictor_sum_squares[i],
            .linear = inputs.rhs[i],
            .scale = inputs.residual_var,
        };

        const auto old_fixed = gelex::detail::compute_posterior_params_core(
            inputs.rhs[i],
            inputs.predictor_sum_squares[i],
            inputs.residual_var,
            fixed_residual_over_prior_var);
        const auto new_fixed = fixed_sampler.posterior_with_logL(kernel);
        Eigen::Vector3d old_fixed_values{
            old_fixed.mean,
            old_fixed.stddev,
            old_fixed.log_likelihood_kernel,
        };
        Eigen::Vector3d new_fixed_values{
            new_fixed.params.mean,
            std::sqrt(new_fixed.params.var),
            new_fixed.log_likelihood_kernel,
        };
        outputs_equivalent
            = outputs_equivalent
              && new_fixed_values.isApprox(old_fixed_values, 1e-13);

        varying_sampler.set_prior_var(inputs.prior_var[i]);
        const auto old_with_log = gelex::detail::compute_posterior_params(
            inputs.rhs[i],
            inputs.prior_var[i],
            inputs.predictor_sum_squares[i],
            inputs.residual_var);
        const auto new_with_log = varying_sampler.posterior_with_logL(kernel);
        Eigen::Vector3d old_with_log_values{
            old_with_log.mean,
            old_with_log.stddev,
            old_with_log.log_likelihood_kernel,
        };
        Eigen::Vector3d new_with_log_values{
            new_with_log.params.mean,
            std::sqrt(new_with_log.params.var),
            new_with_log.log_likelihood_kernel,
        };
        outputs_equivalent
            = outputs_equivalent
              && new_with_log_values.isApprox(old_with_log_values, 1e-13);
    }
    REQUIRE(outputs_equivalent);
}

}  // namespace

TEST_CASE(
    "Normal posterior computation old versus NormalSampler",
    "[!benchmark][stats][normal]")
{
    const auto inputs = make_inputs();
    require_equivalent_outputs(inputs);

    ankerl::nanobench::Bench b;
    b.title("Normal posterior computation")
        .unit("marker")
        .batch(kCount)
        .warmup(5)
        .minEpochIterations(20);

    b.run(
        "old_fixed_prior_core",
        [&]()
        {
            double sink = 0.0;
            const double residual_over_prior_var
                = inputs.residual_var / inputs.fixed_prior_var;
            for (std::size_t i = 0; i < kCount; ++i)
            {
                const auto params
                    = gelex::detail::compute_posterior_params_core(
                        inputs.rhs[i],
                        inputs.predictor_sum_squares[i],
                        inputs.residual_var,
                        residual_over_prior_var);
                sink += params.mean + params.stddev
                        + params.log_likelihood_kernel;
            }
            ankerl::nanobench::doNotOptimizeAway(sink);
        });

    b.run(
        "new_fixed_prior_sampler",
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
        "old_varying_prior",
        [&]()
        {
            double sink = 0.0;
            for (std::size_t i = 0; i < kCount; ++i)
            {
                const auto params = gelex::detail::compute_posterior_params(
                    inputs.rhs[i],
                    inputs.prior_var[i],
                    inputs.predictor_sum_squares[i],
                    inputs.residual_var);
                sink += params.mean + params.stddev
                        + params.log_likelihood_kernel;
            }
            ankerl::nanobench::doNotOptimizeAway(sink);
        });

    b.run(
        "new_varying_prior_sampler",
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
