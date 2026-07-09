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

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <nanobench.h>
#include <random>

namespace
{

constexpr Eigen::Index kMarkers = 1'048'576;

}  // namespace

TEST_CASE(
    "genetic variance multiplier accumulation",
    "[!benchmark][mcmc][genetic][variance]")
{
    Eigen::VectorXd coeffs(kMarkers);
    Eigen::VectorXi assignment(kMarkers);
    Eigen::VectorXd inverse_multiplier(kMarkers);
    Eigen::VectorXd scale(kMarkers);
    Eigen::VectorXd multiplier{{1.0, 0.25, 1.0, 4.0}};
    Eigen::Index active_count = 0;

    std::mt19937_64 rng{42};
    std::normal_distribution<double> normal{0.0, 1.0};
    std::discrete_distribution<int> component{{0.85, 0.05, 0.05, 0.05}};

    for (Eigen::Index i = 0; i < kMarkers; ++i)
    {
        coeffs(i) = normal(rng);
        assignment(i) = component(rng);
        if (assignment(i) != 0)
        {
            ++active_count;
        }
        inverse_multiplier(i)
            = assignment(i) == 0 ? 0.0 : 1.0 / multiplier(assignment(i));
    }

    ankerl::nanobench::Bench b;
    b.title("genetic variance multiplier accumulation")
        .unit("marker")
        .batch(static_cast<double>(kMarkers))
        .warmup(5)
        .minEpochIterations(20);

    b.run(
        "raw_for_loop_baseline",
        [&]
        {
            Eigen::Index n = 0;
            double sum_squares = 0.0;
            for (Eigen::Index i = 0; i < coeffs.size(); ++i)
            {
                const int component_i = assignment(i);
                if (component_i == 0)
                {
                    continue;
                }
                ++n;
                const double coeff = coeffs(i);
                sum_squares += (coeff * coeff) / multiplier(component_i);
            }
            ankerl::nanobench::doNotOptimizeAway(n);
            ankerl::nanobench::doNotOptimizeAway(sum_squares);
        });

    b.run(
        "fill_scale_then_eigen",
        [&]
        {
            Eigen::Index n = 0;
            for (Eigen::Index i = 0; i < coeffs.size(); ++i)
            {
                const int component_i = assignment(i);
                if (component_i == 0)
                {
                    scale(i) = 0.0;
                    continue;
                }
                ++n;
                scale(i) = 1.0 / multiplier(component_i);
            }
            const double sum_squares
                = (coeffs.array().square() * scale.array()).sum();
            ankerl::nanobench::doNotOptimizeAway(n);
            ankerl::nanobench::doNotOptimizeAway(sum_squares);
        });

    b.run(
        "cached_scale_eigen",
        [&]
        {
            const double sum_squares
                = (coeffs.array().square() * inverse_multiplier.array()).sum();
            ankerl::nanobench::doNotOptimizeAway(active_count);
            ankerl::nanobench::doNotOptimizeAway(sum_squares);
        });
}
