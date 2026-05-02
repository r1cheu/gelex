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

#include <cstddef>
#include <vector>

#include <Eigen/Core>

#include <nanobench.h>

#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/infer/detail/marker_op.h"

namespace
{

constexpr Eigen::Index kIndividuals = 50'000;
constexpr std::size_t kMarkers = 256;

struct Inputs
{
    Eigen::VectorXd y_adj;
    Eigen::VectorXd gebv;
    std::vector<Eigen::VectorXd> cols;
    std::vector<double> old_vals;
    std::vector<double> new_vals;
};

auto make_inputs() -> Inputs
{
    Inputs in;
    in.y_adj = Eigen::VectorXd::Zero(kIndividuals);
    in.gebv = Eigen::VectorXd::Zero(kIndividuals);
    in.cols.reserve(kMarkers);
    in.old_vals.resize(kMarkers);
    in.new_vals.resize(kMarkers);

    for (std::size_t marker = 0; marker < kMarkers; ++marker)
    {
        Eigen::VectorXd col(kIndividuals);
        for (Eigen::Index i = 0; i < kIndividuals; ++i)
        {
            const auto v = static_cast<double>((i + (marker * 7)) % 251);
            col(i) = (v * 0.01) - 1.2;
        }
        in.cols.emplace_back(std::move(col));
        in.old_vals[marker] = 0.001 * static_cast<double>((marker % 97) + 1);
        in.new_vals[marker] = 0.001 * static_cast<double>((marker % 53) + 1);
    }

    return in;
}

}  // namespace

TEST_CASE("marker_op apply_marker_update", "[!benchmark][mcmc][marker_op]")
{
    auto inputs = make_inputs();

    ankerl::nanobench::Bench b;
    b.title("apply_marker_update")
        .unit("marker")
        .batch(kMarkers)
        .warmup(5)
        .minEpochIterations(20);

    b.run(
        "swap (old -> new)",
        [&]()
        {
            for (std::size_t marker = 0; marker < kMarkers; ++marker)
            {
                gelex::infer::detail::apply_marker_update(
                    inputs.y_adj,
                    inputs.gebv,
                    {},
                    inputs.cols[marker],
                    {.old_value = inputs.old_vals[marker],
                     .new_value = inputs.new_vals[marker]});
            }
            ankerl::nanobench::doNotOptimizeAway(inputs.y_adj.data());
            ankerl::nanobench::doNotOptimizeAway(inputs.gebv.data());
        });

    b.run(
        "spike (new = 0)",
        [&]()
        {
            for (std::size_t marker = 0; marker < kMarkers; ++marker)
            {
                gelex::infer::detail::apply_marker_update(
                    inputs.y_adj,
                    inputs.gebv,
                    {},
                    inputs.cols[marker],
                    {.old_value = inputs.old_vals[marker], .new_value = 0.0});
            }
            ankerl::nanobench::doNotOptimizeAway(inputs.y_adj.data());
            ankerl::nanobench::doNotOptimizeAway(inputs.gebv.data());
        });
}
