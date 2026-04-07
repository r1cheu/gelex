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

#include <random>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/model/bayes/samplers/detail/stick_breaking.h"
#include "gelex/model/bayes/states.h"

using Catch::Matchers::WithinAbs;
using gelex::bayes::Assignment;
using gelex::detail::sample_stick_posteriors;

TEST_CASE("stick_to_pi produces valid probability vector", "[stick_breaking]")
{
    Eigen::VectorXd q(4);
    q << 0.01, 0.5, 0.8, 0.75;

    auto pi = Assignment::stick_to_pi(q);

    REQUIRE(pi.size() == 5);
    CHECK_THAT(pi.sum(), WithinAbs(1.0, 1e-12));
    for (Eigen::Index i = 0; i < pi.size(); ++i)
    {
        CHECK(pi(i) >= 0.0);
    }
}

TEST_CASE("stick_to_pi K=1 (SpikePrior)", "[stick_breaking]")
{
    Eigen::VectorXd q(1);
    q << 0.05;

    auto pi = Assignment::stick_to_pi(q);

    REQUIRE(pi.size() == 2);
    CHECK_THAT(pi(0), WithinAbs(0.95, 1e-12));
    CHECK_THAT(pi(1), WithinAbs(0.05, 1e-12));
}

TEST_CASE("round-trip pi_to_stick(stick_to_pi(q))", "[stick_breaking]")
{
    Eigen::VectorXd q(4);
    q << 0.01, 0.5, 0.8, 0.75;

    auto recovered = Assignment::pi_to_stick(Assignment::stick_to_pi(q));

    REQUIRE(recovered.size() == q.size());
    CHECK(q.isApprox(recovered, 1e-12));
}

TEST_CASE("round-trip stick_to_pi(pi_to_stick(pi))", "[stick_breaking]")
{
    Eigen::VectorXd pi(5);
    pi << 0.99, 0.005, 0.002, 0.002, 0.001;

    auto recovered = Assignment::stick_to_pi(Assignment::pi_to_stick(pi));

    REQUIRE(recovered.size() == pi.size());
    CHECK(pi.isApprox(recovered, 1e-12));
}

TEST_CASE("sample_stick_posteriors produces valid q", "[stick_breaking]")
{
    std::mt19937_64 rng{42};
    Eigen::VectorXi counts(5);
    counts << 900, 50, 30, 15, 5;

    auto q = sample_stick_posteriors(counts, rng);

    REQUIRE(q.size() == 4);
    for (Eigen::Index k = 0; k < q.size(); ++k)
    {
        CHECK(q(k) > 0.0);
        CHECK(q(k) < 1.0);
    }

    auto pi = Assignment::stick_to_pi(q);
    CHECK_THAT(pi.sum(), WithinAbs(1.0, 1e-12));
}

TEST_CASE("sample_stick_posteriors K=1", "[stick_breaking]")
{
    std::mt19937_64 rng{42};
    Eigen::VectorXi counts(2);
    counts << 950, 50;

    auto q = sample_stick_posteriors(counts, rng);

    REQUIRE(q.size() == 1);
    CHECK(q(0) > 0.0);
    CHECK(q(0) < 1.0);
}

TEST_CASE(
    "Assignment initializes stick_probs from proportion",
    "[stick_breaking]")
{
    Eigen::VectorXd init_pi(5);
    init_pi << 0.99, 0.005, 0.002, 0.002, 0.001;

    Assignment asgn(10, init_pi);

    REQUIRE(asgn.stick_probs.size() == 4);
    auto recovered_pi = Assignment::stick_to_pi(asgn.stick_probs);
    CHECK(init_pi.isApprox(recovered_pi, 1e-12));
}
