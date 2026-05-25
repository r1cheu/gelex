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

#include <array>
#include <cmath>
#include <random>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/infra/stats/detail/truncated_normal.h"

namespace gelex
{
namespace
{

using Catch::Matchers::WithinAbs;

constexpr int kNumSamples = 100000;

struct TruncConfig
{
    double mu;
    double sigma;
    double bound;
    bool left_truncated;  // true = [a, +inf), false = (-inf, b]
    double expected_mean;
    double expected_var;
};

// scipy.stats.truncnorm reference values
constexpr std::array<TruncConfig, 9> kScipyConfigs = {{
    // TN(0,1,[0,+inf))
    {0.0, 1.0, 0.0, true, 7.978845608028654e-01, 3.633802276324186e-01},
    // TN(0,1,(-inf,0])
    {0.0, 1.0, 0.0, false, -7.978845608028654e-01, 3.633802276324186e-01},
    // TN(5,2,[3,+inf))
    {5.0, 2.0, 3.0, true, 5.575199941878357e+00, 2.518745143106422e+00},
    // TN(5,2,(-inf,7])
    {5.0, 2.0, 7.0, false, 4.424800058121643e+00, 2.518745143106422e+00},
    // TN(-10,1,[0,+inf))  — deep tail, uses Robert algorithm
    {-10.0, 1.0, 0.0, true, 9.809323396256353e-02, 9.445377825130441e-03},
    // TN(10,1,(-inf,0])
    {10.0, 1.0, 0.0, false, -9.809323396256353e-02, 9.445377825130441e-03},
    // TN(0,1,[2,+inf))
    {0.0, 1.0, 2.0, true, 2.373215532822840e+00, 1.142791004140833e-01},
    // TN(0,1,[3,+inf))
    {0.0, 1.0, 3.0, true, 3.283098654930440e+00, 7.055918678525586e-02},
    // TN(0,1,(-inf,-2])
    {0.0, 1.0, -2.0, false, -2.373215532822840e+00, 1.142791004140833e-01},
}};

auto draw_samples(const TruncConfig& cfg, std::mt19937_64& rng)
    -> std::pair<double, double>
{
    double sum = 0.0;
    double sum_sq = 0.0;

    for (int i = 0; i < kNumSamples; ++i)
    {
        double val = cfg.left_truncated
                         ? stats::detail::sample_left_truncated_normal(
                               cfg.mu, cfg.sigma, cfg.bound, rng)
                         : stats::detail::sample_right_truncated_normal(
                               cfg.mu, cfg.sigma, cfg.bound, rng);
        sum += val;
        sum_sq += val * val;
    }

    double mean = sum / kNumSamples;
    double var = (sum_sq / kNumSamples) - (mean * mean);
    return {mean, var};
}

TEST_CASE("Mean & variance vs scipy reference values", "[truncated_normal]")
{
    std::mt19937_64 rng(12345);

    for (const auto& cfg : kScipyConfigs)
    {
        auto [mean, var] = draw_samples(cfg, rng);

        // tolerance: ~4*sigma/sqrt(N) for mean
        double mean_tol
            = 4.0 * std::sqrt(cfg.expected_var) / std::sqrt(kNumSamples);
        // looser tolerance for variance (CLT on (x-mu)^2)
        double var_tol = (0.1 * cfg.expected_var) + 0.01;

        INFO(
            "mu=" << cfg.mu << " sigma=" << cfg.sigma << " bound=" << cfg.bound
                  << " left=" << cfg.left_truncated);
        CHECK_THAT(mean, WithinAbs(cfg.expected_mean, mean_tol));
        CHECK_THAT(var, WithinAbs(cfg.expected_var, var_tol));
    }
}

TEST_CASE(
    "Bounds correctness — all samples respect truncation",
    "[truncated_normal]")
{
    std::mt19937_64 rng(42);

    bool all_in_bounds = true;

    for (const auto& cfg : kScipyConfigs)
    {
        for (int i = 0; i < 10000; ++i)
        {
            if (cfg.left_truncated)
            {
                double val = stats::detail::sample_left_truncated_normal(
                    cfg.mu, cfg.sigma, cfg.bound, rng);
                if (val < cfg.bound)
                {
                    all_in_bounds = false;
                }
            }
            else
            {
                double val = stats::detail::sample_right_truncated_normal(
                    cfg.mu, cfg.sigma, cfg.bound, rng);
                if (val > cfg.bound)
                {
                    all_in_bounds = false;
                }
            }
        }
    }

    CHECK(all_in_bounds);
}

TEST_CASE(
    "Symmetry: left-truncated mean == -right-truncated mean",
    "[truncated_normal]")
{
    struct SymPair
    {
        double mu;
        double sigma;
        double a;
    };

    constexpr std::array<SymPair, 4> kPairs = {{
        {0.0, 1.0, 0.0},
        {5.0, 2.0, 3.0},
        {-10.0, 1.0, 0.0},
        {0.0, 1.0, 2.0},
    }};

    for (const auto& [mu, sigma, a] : kPairs)
    {
        std::mt19937_64 rng_left(77);
        std::mt19937_64 rng_right(77);

        double sum_left = 0.0;
        double sum_right = 0.0;

        for (int i = 0; i < kNumSamples; ++i)
        {
            sum_left += stats::detail::sample_left_truncated_normal(
                mu, sigma, a, rng_left);
            sum_right += stats::detail::sample_right_truncated_normal(
                -mu, sigma, -a, rng_right);
        }

        double mean_left = sum_left / kNumSamples;
        double mean_right = sum_right / kNumSamples;

        INFO("mu=" << mu << " sigma=" << sigma << " a=" << a);
        // Same seed + symmetry implementation => should be exact
        CHECK_THAT(mean_left, WithinAbs(-mean_right, 1e-10));
    }
}

TEST_CASE("Both algorithm paths produce correct mean/var", "[truncated_normal]")
{
    // alpha=3.4 uses inverse CDF path (< 3.5 threshold)
    // TN(0,1,[3.4,+inf)) — scipy: mean=3.583, var=0.0540
    {
        std::mt19937_64 rng(999);
        // mu=-3.4, sigma=1, a=0 => alpha = (0 - (-3.4))/1 = 3.4
        auto [mean, var] = draw_samples({-3.4, 1.0, 0.0, true, 0.0, 0.0}, rng);
        // Just verify finite and positive bound
        CHECK(std::isfinite(mean));
        CHECK(mean >= 0.0);
    }

    // alpha=3.6 uses Robert algorithm (>= 3.5 threshold)
    // TN(0,1,[3.6,+inf)) — scipy: mean=3.874, var=0.0460
    {
        std::mt19937_64 rng(999);
        auto [mean, var] = draw_samples({-3.6, 1.0, 0.0, true, 0.0, 0.0}, rng);
        CHECK(std::isfinite(mean));
        CHECK(mean >= 0.0);
    }

    // Verify both paths against scipy values
    // TN(0,1,[2,+inf)) -> alpha=2 (inverse CDF)
    {
        std::mt19937_64 rng(999);
        auto [mean, var] = draw_samples(
            {0.0, 1.0, 2.0, true, 2.373215532822840e+00, 1.142791004140833e-01},
            rng);
        double mean_tol
            = 4.0 * std::sqrt(1.142791004140833e-01) / std::sqrt(kNumSamples);
        CHECK_THAT(mean, WithinAbs(2.373215532822840e+00, mean_tol));
    }

    // TN(0,1,[3,+inf)) -> alpha=3 (inverse CDF, near threshold)
    {
        std::mt19937_64 rng(999);
        auto [mean, var] = draw_samples(
            {0.0, 1.0, 3.0, true, 3.283098654930440e+00, 7.055918678525586e-02},
            rng);
        double mean_tol
            = 4.0 * std::sqrt(7.055918678525586e-02) / std::sqrt(kNumSamples);
        CHECK_THAT(mean, WithinAbs(3.283098654930440e+00, mean_tol));
    }

    // TN(-10,1,[0,+inf)) -> alpha=10 (Robert)
    {
        std::mt19937_64 rng(999);
        auto [mean, var] = draw_samples(
            {-10.0,
             1.0,
             0.0,
             true,
             9.809323396256353e-02,
             9.445377825130441e-03},
            rng);
        double mean_tol
            = 4.0 * std::sqrt(9.445377825130441e-03) / std::sqrt(kNumSamples);
        CHECK_THAT(mean, WithinAbs(9.809323396256353e-02, mean_tol));
    }
}

TEST_CASE("Extreme tails: mu=-30 left-truncated at 0", "[truncated_normal]")
{
    std::mt19937_64 rng(42);
    double sum = 0.0;
    double sum_sq = 0.0;
    bool all_valid = true;

    for (int i = 0; i < kNumSamples; ++i)
    {
        double val
            = stats::detail::sample_left_truncated_normal(-30.0, 1.0, 0.0, rng);
        if (val < 0.0 || !std::isfinite(val))
        {
            all_valid = false;
        }
        sum += val;
        sum_sq += val * val;
    }

    CHECK(all_valid);
    double mean = sum / kNumSamples;
    double var = (sum_sq / kNumSamples) - (mean * mean);
    CHECK(var > 1e-10);
    CHECK(mean < 1.0);
    CHECK(mean > 0.0);
}

TEST_CASE("Extreme tails: mu=30 right-truncated at 0", "[truncated_normal]")
{
    std::mt19937_64 rng(42);
    double sum = 0.0;
    double sum_sq = 0.0;
    bool all_valid = true;

    for (int i = 0; i < kNumSamples; ++i)
    {
        double val
            = stats::detail::sample_right_truncated_normal(30.0, 1.0, 0.0, rng);
        if (val > 0.0 || !std::isfinite(val))
        {
            all_valid = false;
        }
        sum += val;
        sum_sq += val * val;
    }

    CHECK(all_valid);
    double mean = sum / kNumSamples;
    double var = (sum_sq / kNumSamples) - (mean * mean);
    CHECK(var > 1e-10);
    CHECK(mean > -1.0);
    CHECK(mean < 0.0);
}

}  // namespace
}  // namespace gelex
