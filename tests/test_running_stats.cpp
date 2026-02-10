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
#include <limits>

#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>

#include "gelex/exception.h"
#include "gelex/utils/running_stats.h"

namespace gelex
{
namespace
{

constexpr double k_tolerance = 1e-12;

auto sample_stddev(Eigen::Ref<const Eigen::VectorXd> values) -> double
{
    if (values.size() <= 1)
    {
        return 0.0;
    }

    double mean = values.mean();
    double sum_sq_diff = (values.array() - mean).square().sum();
    return std::sqrt(sum_sq_diff / static_cast<double>(values.size() - 1));
}

auto require_vector_is_approx(
    Eigen::Ref<const Eigen::VectorXd> actual,
    Eigen::Ref<const Eigen::VectorXd> expected,
    double tolerance = k_tolerance) -> void
{
    REQUIRE(actual.size() == expected.size());
    REQUIRE(actual.isApprox(expected, tolerance));
}

auto compute_row_mean(Eigen::Ref<const Eigen::MatrixXd> matrix)
    -> Eigen::VectorXd
{
    Eigen::VectorXd mean(matrix.rows());
    for (Eigen::Index i = 0; i < matrix.rows(); ++i)
    {
        mean(i) = matrix.row(i).mean();
    }
    return mean;
}

auto compute_row_sample_stddev(Eigen::Ref<const Eigen::MatrixXd> matrix)
    -> Eigen::VectorXd
{
    Eigen::VectorXd stddev(matrix.rows());
    for (Eigen::Index i = 0; i < matrix.rows(); ++i)
    {
        stddev(i) = sample_stddev(matrix.row(i).transpose());
    }
    return stddev;
}

}  // namespace

TEST_CASE("RunningStats default state is empty", "[utils][running_stats]")
{
    RunningStats stats;

    RunningStatsResult result = stats.result();
    REQUIRE(result.mean.size() == 0);
    REQUIRE(result.stddev.size() == 0);
}

TEST_CASE(
    "RunningStats computes row-wise mean and stddev",
    "[utils][running_stats]")
{
    RunningStats stats;
    Eigen::MatrixXd block(2, 4);
    block << 1.0, 2.0, 3.0, 4.0, 10.0, 20.0, 30.0, 40.0;

    stats.update(block);

    RunningStatsResult result = stats.result();
    Eigen::VectorXd expected_mean(2);
    expected_mean << 2.5, 25.0;
    Eigen::VectorXd expected_stddev(2);
    expected_stddev << std::sqrt(5.0 / 3.0), std::sqrt(500.0 / 3.0);

    require_vector_is_approx(result.mean, expected_mean);
    require_vector_is_approx(result.stddev, expected_stddev);
}

TEST_CASE(
    "RunningStats count=1 gives zero row stddev",
    "[utils][running_stats]")
{
    RunningStats stats;
    Eigen::MatrixXd block(3, 1);
    block << 3.0, -1.5, 0.25;

    stats.update(block);

    RunningStatsResult result = stats.result();
    Eigen::VectorXd expected_mean(3);
    expected_mean << 3.0, -1.5, 0.25;
    Eigen::VectorXd expected_stddev = Eigen::VectorXd::Zero(3);

    require_vector_is_approx(result.mean, expected_mean);
    require_vector_is_approx(result.stddev, expected_stddev);
}

TEST_CASE(
    "RunningStats batched updates match one-shot for axis=1",
    "[utils][running_stats]")
{
    RunningStats batched;
    RunningStats one_shot;

    Eigen::MatrixXd part1(3, 2);
    part1 << 1.0, 2.0, 4.0, 5.0, 7.0, 8.0;
    Eigen::MatrixXd part2(3, 3);
    part2 << 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0;
    Eigen::MatrixXd all(3, 5);
    all << 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 5.0, 6.0, 7.0, 8.0, 7.0, 8.0, 9.0,
        10.0, 11.0;

    batched.update(part1);
    batched.update(part2);
    one_shot.update(all);

    RunningStatsResult batched_result = batched.result();
    RunningStatsResult one_shot_result = one_shot.result();
    require_vector_is_approx(batched_result.mean, one_shot_result.mean);
    require_vector_is_approx(batched_result.stddev, one_shot_result.stddev);
}

TEST_CASE("RunningStats empty-column update is no-op", "[utils][running_stats]")
{
    RunningStats stats;
    Eigen::MatrixXd block(2, 3);
    block << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    Eigen::MatrixXd empty_block(2, 0);

    stats.update(block);
    RunningStatsResult before = stats.result();

    stats.update(empty_block);
    RunningStatsResult after = stats.result();

    require_vector_is_approx(after.mean, before.mean);
    require_vector_is_approx(after.stddev, before.stddev);
}

TEST_CASE(
    "RunningStats rejects row-size mismatch across updates",
    "[utils][running_stats]")
{
    RunningStats stats;
    Eigen::MatrixXd valid(2, 2);
    valid << 1.0, 2.0, 3.0, 4.0;
    Eigen::MatrixXd mismatch(3, 2);
    mismatch << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

    stats.update(valid);
    RunningStatsResult before = stats.result();

    REQUIRE_THROWS_AS(stats.update(mismatch), InvalidInputException);

    RunningStatsResult after = stats.result();
    require_vector_is_approx(after.mean, before.mean);
    require_vector_is_approx(after.stddev, before.stddev);
}

TEST_CASE("RunningStats rejects NaN and Inf", "[utils][running_stats]")
{
    RunningStats stats;

    Eigen::MatrixXd with_nan(2, 2);
    with_nan << 1.0, std::numeric_limits<double>::quiet_NaN(), 3.0, 4.0;
    Eigen::MatrixXd with_inf(2, 2);
    with_inf << 1.0, std::numeric_limits<double>::infinity(), 3.0, 4.0;
    Eigen::MatrixXd with_neg_inf(2, 2);
    with_neg_inf << 1.0, -std::numeric_limits<double>::infinity(), 3.0, 4.0;

    REQUIRE_THROWS_AS(stats.update(with_nan), InvalidInputException);
    REQUIRE_THROWS_AS(stats.update(with_inf), InvalidInputException);
    REQUIRE_THROWS_AS(stats.update(with_neg_inf), InvalidInputException);
}

TEST_CASE(
    "RunningStats invalid first update keeps empty state",
    "[utils][running_stats]")
{
    RunningStats stats;

    Eigen::MatrixXd with_nan(2, 2);
    with_nan << 1.0, std::numeric_limits<double>::quiet_NaN(), 3.0, 4.0;

    REQUIRE_THROWS_AS(stats.update(with_nan), InvalidInputException);

    RunningStatsResult result = stats.result();
    REQUIRE(result.mean.size() == 0);
    REQUIRE(result.stddev.size() == 0);
}

TEST_CASE(
    "RunningStats rejects zero-row non-empty updates",
    "[utils][running_stats]")
{
    RunningStats stats;
    Eigen::MatrixXd zero_rows(0, 2);

    REQUIRE_THROWS_AS(stats.update(zero_rows), InvalidInputException);

    RunningStatsResult result = stats.result();
    REQUIRE(result.mean.size() == 0);
    REQUIRE(result.stddev.size() == 0);
}

TEST_CASE(
    "RunningStats exception keeps state unchanged",
    "[utils][running_stats]")
{
    RunningStats stats;
    Eigen::MatrixXd valid(2, 3);
    valid << 1.0, 2.0, 3.0, 7.0, 8.0, 9.0;

    stats.update(valid);
    RunningStatsResult before = stats.result();

    Eigen::MatrixXd invalid(2, 2);
    invalid << 4.0, std::numeric_limits<double>::quiet_NaN(), 10.0, 11.0;
    REQUIRE_THROWS_AS(stats.update(invalid), InvalidInputException);

    RunningStatsResult after = stats.result();
    require_vector_is_approx(after.mean, before.mean);
    require_vector_is_approx(after.stddev, before.stddev);
}

TEST_CASE(
    "RunningStats supports mixed input element types",
    "[utils][running_stats]")
{
    RunningStats stats;

    Eigen::MatrixXf float_block(2, 2);
    float_block << 1.0f, 2.0f, 10.0f, 20.0f;
    Eigen::MatrixXi int_block(2, 2);
    int_block << 3, 4, 30, 40;
    Eigen::MatrixXd double_block(2, 2);
    double_block << 5.0, 6.0, 50.0, 60.0;

    stats.update(float_block);
    stats.update(int_block);
    stats.update(double_block);

    Eigen::MatrixXd full(2, 6);
    full << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0;
    RunningStatsResult result = stats.result();

    require_vector_is_approx(result.mean, compute_row_mean(full), 1e-10);
    require_vector_is_approx(
        result.stddev, compute_row_sample_stddev(full), 1e-10);
}

TEST_CASE(
    "RunningStats is numerically stable for large row values",
    "[utils][running_stats]")
{
    RunningStats stats;
    Eigen::MatrixXd block(2, 3);
    block << 1.0e12, 1.0e12 + 1.0, 1.0e12 + 2.0, 5.0e11, 5.0e11 + 2.0,
        5.0e11 + 4.0;

    stats.update(block);

    RunningStatsResult result = stats.result();
    REQUIRE(result.mean.size() == 2);
    REQUIRE(result.stddev.size() == 2);
    for (Eigen::Index i = 0; i < result.mean.size(); ++i)
    {
        REQUIRE(std::isfinite(result.mean(i)));
        REQUIRE(std::isfinite(result.stddev(i)));
    }

    require_vector_is_approx(result.mean, compute_row_mean(block), 1e-6);
    require_vector_is_approx(
        result.stddev, compute_row_sample_stddev(block), 1e-9);
}

}  // namespace gelex
