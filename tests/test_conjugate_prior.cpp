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

#include <cstdint>
#include <random>

#include <Eigen/Core>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/exception.h"
#include "gelex/infra/stats/conjugate_prior.h"

namespace gelex
{
namespace
{

constexpr std::uint64_t kSeed = 0xC0FFEE1234ULL;
constexpr int kDrawCount = 20000;

TEST_CASE("BetaSampler validates input", "[conjugate][beta]")
{
    REQUIRE_THROWS_AS(BetaSampler<double>(0.0, 1.0), GelexException);
    REQUIRE_THROWS_AS(BetaSampler<double>(1.0, 0.0), GelexException);
    REQUIRE_THROWS_AS(BetaSampler<double>(-1.0, 1.0), GelexException);

    BetaSampler<double> s{1.0, 1.0};
    std::mt19937_64 rng{kSeed};
    REQUIRE_THROWS_AS(s({-1, 0}, rng), GelexException);
    REQUIRE_THROWS_AS(s({0, -1}, rng), GelexException);
}

TEMPLATE_TEST_CASE(
    "BetaSampler posterior mean concentrates near MLE for large samples",
    "[conjugate][beta]",
    float,
    double)
{
    using T = TestType;
    BetaSampler<T> s{T{1}, T{1}};
    typename BetaSampler<T>::Likelihood lik{800, 200};
    std::mt19937_64 rng{kSeed};

    Eigen::VectorX<T> draws(kDrawCount);
    for (int i = 0; i < kDrawCount; ++i)
    {
        draws(i) = s(lik, rng);
    }

    REQUIRE((draws.array() > T{0}).all());
    REQUIRE((draws.array() < T{1}).all());

    const T mean = draws.mean();
    const T var
        = ((draws.array() - mean).square()).sum() / static_cast<T>(kDrawCount);

    REQUIRE_THAT(
        static_cast<double>(mean),
        Catch::Matchers::WithinAbs(801.0 / 1002.0, 5e-3));
    // Beta(α=801, β=201) variance = αβ / ((α+β)² (α+β+1))
    REQUIRE_THAT(
        static_cast<double>(var), Catch::Matchers::WithinAbs(1.5988e-4, 1e-5));
}

TEST_CASE("BetaSampler is reproducible under fixed seed", "[conjugate][beta]")
{
    BetaSampler<double> s{2.0, 3.0};
    BetaSampler<double>::Likelihood lik{10, 5};
    std::mt19937_64 r1{kSeed};
    std::mt19937_64 r2{kSeed};

    Eigen::VectorXd a(100);
    Eigen::VectorXd b(100);
    for (int i = 0; i < 100; ++i)
    {
        a(i) = s(lik, r1);
        b(i) = s(lik, r2);
    }
    REQUIRE(a.isApprox(b));
}

TEST_CASE("DirichletSampler validates input", "[conjugate][dirichlet]")
{
    REQUIRE_THROWS_AS(
        DirichletSampler<double>(Eigen::VectorXd{{1.0}}), GelexException);
    REQUIRE_THROWS_AS(
        DirichletSampler<double>(Eigen::VectorXd{{1.0, 0.0, 1.0}}),
        GelexException);

    DirichletSampler<double> s{Eigen::VectorXd{{1.0, 1.0, 1.0}}};

    std::mt19937_64 rng{kSeed};
    REQUIRE_THROWS_AS(s(Eigen::VectorXi{{1, 2}}, rng), GelexException);
    REQUIRE_THROWS_AS(s(Eigen::VectorXi{{-1, 2, 3}}, rng), GelexException);
}

TEMPLATE_TEST_CASE(
    "DirichletSampler posterior mean and simplex invariants",
    "[conjugate][dirichlet]",
    float,
    double)
{
    using T = TestType;
    const Eigen::VectorX<T> alpha{{T{1}, T{1}, T{1}}};
    DirichletSampler<T> s{alpha};

    const Eigen::VectorXi counts{{600, 300, 100}};
    std::mt19937_64 rng{kSeed};

    Eigen::Matrix<T, 3, Eigen::Dynamic> draws(3, kDrawCount);
    for (int i = 0; i < kDrawCount; ++i)
    {
        draws.col(i) = s(counts, rng);
    }

    // simplex invariants
    REQUIRE((draws.array() >= T{0}).all());
    REQUIRE((draws.array() <= T{1}).all());
    const Eigen::RowVectorX<T> sums = draws.colwise().sum();
    const Eigen::RowVectorX<T> ones = Eigen::RowVectorX<T>::Ones(kDrawCount);
    REQUIRE(sums.isApprox(ones));

    // posterior mean matches Dirichlet(α + counts) mean
    const Eigen::VectorX<T> mean = draws.rowwise().mean();
    const Eigen::VectorX<T> posterior_alpha = alpha + counts.cast<T>();
    const Eigen::VectorX<T> expected = posterior_alpha / posterior_alpha.sum();
    REQUIRE(mean.isApprox(expected, T{5e-3}));
}

TEST_CASE(
    "DirichletSampler is reproducible under fixed seed",
    "[conjugate][dirichlet]")
{
    constexpr int kIters = 50;
    DirichletSampler<double> s{Eigen::VectorXd{{1.5, 2.0, 0.5}}};
    const Eigen::VectorXi counts{{4, 7, 2}};

    std::mt19937_64 r1{kSeed};
    std::mt19937_64 r2{kSeed};
    Eigen::MatrixXd a(3, kIters);
    Eigen::MatrixXd b(3, kIters);
    for (int i = 0; i < kIters; ++i)
    {
        a.col(i) = s(counts, r1);
        b.col(i) = s(counts, r2);
    }
    REQUIRE(a.isApprox(b));
}

TEST_CASE("ScaledInvChi2Sampler validates input", "[conjugate][invchi2]")
{
    REQUIRE_THROWS_AS(ScaledInvChi2Sampler<double>(1.0, -1.0), GelexException);

    ScaledInvChi2Sampler<double> s{2.0, 1.0};
    std::mt19937_64 rng{kSeed};
    REQUIRE_THROWS_AS(s({-1, 0.0}, rng), GelexException);
    REQUIRE_THROWS_AS(s({1, -1.0}, rng), GelexException);

    // Reference improper prior (nu0=-2, s2_0=0): posterior nu1 = n - 2 must
    // be > 0, so n <= 2 is rejected.
    ScaledInvChi2Sampler<double> ref{-2.0, 0.0};
    REQUIRE_THROWS_AS(ref({2, 5.0}, rng), GelexException);

    // Negative posterior s2: nu0=-1, s2_0=10 → nu0*s2_0 = -10, with
    // sum_squares=0 the posterior s2 is negative.
    ScaledInvChi2Sampler<double> neg{-1.0, 10.0};
    REQUIRE_THROWS_AS(neg({2, 0.0}, rng), GelexException);
}

TEMPLATE_TEST_CASE(
    "ScaledInvChi2Sampler posterior mean approaches σ²_true with much data",
    "[conjugate][invchi2]",
    float,
    double)
{
    using T = TestType;
    ScaledInvChi2Sampler<T> s{T{4}, T{1}};
    constexpr Eigen::Index kN = 5000;
    constexpr T kTrueVar = T{4};
    typename ScaledInvChi2Sampler<T>::Likelihood lik{
        kN, kTrueVar * static_cast<T>(kN)};
    std::mt19937_64 rng{kSeed};

    Eigen::VectorX<T> draws(kDrawCount);
    for (int i = 0; i < kDrawCount; ++i)
    {
        draws(i) = s(lik, rng);
    }

    REQUIRE((draws.array() > T{0}).all());

    const T mean = draws.mean();
    const T nu1 = T{4} + static_cast<T>(kN);
    const T s2_1 = ((T{4} * T{1}) + (kTrueVar * static_cast<T>(kN))) / nu1;
    const T expected = (nu1 * s2_1) / (nu1 - T{2});
    REQUIRE_THAT(
        static_cast<double>(mean),
        Catch::Matchers::WithinAbs(static_cast<double>(expected), 5e-2));
}

TEST_CASE(
    "ScaledInvChi2Sampler is reproducible under fixed seed",
    "[conjugate][invchi2]")
{
    constexpr int kIters = 100;
    ScaledInvChi2Sampler<double> s{3.0, 0.5};
    ScaledInvChi2Sampler<double>::Likelihood lik{20, 12.5};

    std::mt19937_64 r1{kSeed};
    std::mt19937_64 r2{kSeed};
    Eigen::VectorXd a(kIters);
    Eigen::VectorXd b(kIters);
    for (int i = 0; i < kIters; ++i)
    {
        a(i) = s(lik, r1);
        b(i) = s(lik, r2);
    }
    REQUIRE(a.isApprox(b));
}

}  // namespace
}  // namespace gelex
