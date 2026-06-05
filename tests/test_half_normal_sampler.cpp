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

#include <algorithm>
#include <catch2/matchers/catch_matchers.hpp>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <random>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/infra/stats/half_normal_sampler.h"
#include "gelex/infra/stats/normal_sampler.h"

namespace gelex
{

using stats::HalfNormalSampler;
using stats::NormalSampler;

namespace
{

constexpr std::uint64_t SEED = 0xDEADBEEF42ULL;
constexpr int DRAW_COUNT = 1'000'000;

TEST_CASE(
    "HalfNormalSampler::posterior matches NormalSampler::posterior",
    "[stats][half_normal_sampler]")
{
    constexpr double PRIOR_VAR = 2.5;
    HalfNormalSampler<double> hn{PRIOR_VAR};
    NormalSampler<double> n{PRIOR_VAR};

    const NormalSampler<double>::Kernel kernel{
        .quadratic = 3.0,
        .linear = 1.5,
        .scale = 1.0,
    };

    const auto hn_params = hn.posterior(kernel);
    const auto n_params = n.posterior(kernel);

    REQUIRE_THAT(
        hn_params.mean, Catch::Matchers::WithinRel(n_params.mean, 1e-12));
    REQUIRE_THAT(
        hn_params.var, Catch::Matchers::WithinRel(n_params.var, 1e-12));
}

TEST_CASE(
    "HalfNormalSampler log_marginal_kernel includes sign tail mass",
    "[stats][half_normal_sampler]")
{
    constexpr double PRIOR_VAR = 3.0;
    HalfNormalSampler<double> hn{PRIOR_VAR};
    NormalSampler<double> n{PRIOR_VAR};

    const NormalSampler<double>::Kernel kernel{
        .quadratic = 4.0,
        .linear = 2.0,
        .scale = 1.5,
    };

    const auto hn_pos = hn.posterior_with_logL(kernel, +1);
    const auto hn_neg = hn.posterior_with_logL(kernel, -1);
    const auto n_post = n.posterior_with_logL(kernel);
    const double z = n_post.params.mean / std::sqrt(n_post.params.var);
    const double log_tail_pos
        = std::log(0.5) + std::log(std::erfc(-z / std::numbers::sqrt2));
    const double log_tail_neg
        = std::log(0.5) + std::log(std::erfc(z / std::numbers::sqrt2));
    const double normal_half_log
        = n_post.log_likelihood_kernel + std::log(2.0);

    REQUIRE_THAT(
        hn_pos.log_marginal_kernel,
        Catch::Matchers::WithinAbs(normal_half_log + log_tail_pos, 1e-12));
    REQUIRE(hn_pos.sign == +1);

    REQUIRE_THAT(
        hn_neg.log_marginal_kernel,
        Catch::Matchers::WithinAbs(normal_half_log + log_tail_neg, 1e-12));
    REQUIRE(hn_neg.sign == -1);

    const double big
        = std::max(hn_pos.log_marginal_kernel, hn_neg.log_marginal_kernel);
    const double logsumexp
        = big
          + std::log(
              std::exp(hn_pos.log_marginal_kernel - big)
              + std::exp(hn_neg.log_marginal_kernel - big));
    REQUIRE_THAT(logsumexp, Catch::Matchers::WithinAbs(normal_half_log, 1e-12));
}

TEST_CASE(
    "HalfNormalSampler log_tail values are finite for large |z| and match "
    "asymptotic",
    "[stats][half_normal_sampler]")
{
    constexpr double SQRT_2 = std::numbers::sqrt2;

    for (const double sign_val : {1.0, -1.0})
    {
        HalfNormalSampler<double> hn{1.0};
        NormalSampler<double> n{1.0};
        const NormalSampler<double>::Kernel kernel{
            .quadratic = 1.0,
            .linear = sign_val * 10.0 * SQRT_2,
            .scale = 1.0,
        };
        const auto pos = hn.posterior_with_logL(kernel, +1);
        const auto neg = hn.posterior_with_logL(kernel, -1);
        const auto normal_post = n.posterior_with_logL(kernel);
        const double z = normal_post.params.mean / std::sqrt(normal_post.params.var);
        const double normal_half_log
            = normal_post.log_likelihood_kernel + std::log(2.0);
        const double log_tail_pos = pos.log_marginal_kernel - normal_half_log;
        const double log_tail_neg = neg.log_marginal_kernel - normal_half_log;

        REQUIRE(std::isfinite(pos.log_marginal_kernel));
        REQUIRE(std::isfinite(neg.log_marginal_kernel));
        REQUIRE(!std::isnan(pos.log_marginal_kernel));
        REQUIRE(!std::isnan(neg.log_marginal_kernel));

        auto log_phi_asymp = [](double zv) -> double
        {
            const double z2 = zv * zv;
            const double z4 = z2 * z2;
            const double z6 = z4 * z2;
            return (-0.5 * z2) - (0.5 * std::log(2.0 * std::numbers::pi))
                   - std::log(std::abs(zv))
                   + std::log1p((-1.0 / z2) + (3.0 / z4) - (15.0 / z6));
        };

        if (sign_val < 0)
        {
            const double ref = log_phi_asymp(z);
            REQUIRE_THAT(
                log_tail_pos, Catch::Matchers::WithinAbs(ref, 1e-6));
        }
        else
        {
            // z ≈ +10: log_tail_neg ≈ log_phi_asymp(-10)
            const double ref = log_phi_asymp(-z);
            REQUIRE_THAT(
                log_tail_neg, Catch::Matchers::WithinAbs(ref, 1e-6));
        }
    }
}

TEST_CASE(
    "HalfNormalSampler signed log marginal sums to sign-free half-normal mass",
    "[stats][half_normal_sampler]")
{
    HalfNormalSampler<double> hn{1.0};
    NormalSampler<double> n{1.0};

    for (const double mean_val : {-5.0, -1.0, 0.0, 1.0, 5.0})
    {
        const NormalSampler<double>::Kernel kernel{
            .quadratic = 1.0,
            .linear = mean_val * 2.0,
            .scale = 1.0,
        };
        const auto pos = hn.posterior_with_logL(kernel, +1);
        const auto neg = hn.posterior_with_logL(kernel, -1);
        const auto normal_post = n.posterior_with_logL(kernel);

        const double lp = pos.log_marginal_kernel;
        const double ln = neg.log_marginal_kernel;
        const double big = std::max(lp, ln);
        const double logsumexp
            = big + std::log(std::exp(lp - big) + std::exp(ln - big));

        REQUIRE_THAT(
            logsumexp,
            Catch::Matchers::WithinAbs(
                normal_post.log_likelihood_kernel + std::log(2.0), 1e-12));
    }
}

TEST_CASE(
    "HalfNormalSampler draw half-normal moments sign=+1 mean=0 var=1",
    "[stats][half_normal_sampler]")
{
    HalfNormalSampler<double> hn{1.0};
    const HalfNormalSampler<double>::Params post{.mean = 0.0, .var = 1.0};

    std::mt19937_64 rng{SEED};
    double sum = 0.0;
    double sum2 = 0.0;
    for (int i = 0; i < DRAW_COUNT; ++i)
    {
        const double x = hn.draw(post, +1, rng);
        sum += x;
        sum2 += x * x;
    }
    const double mean = sum / DRAW_COUNT;
    const double var = (sum2 / DRAW_COUNT) - (mean * mean);

    const double expected_mean = std::sqrt(2.0 / std::numbers::pi);
    const double expected_var = 1.0 - (2.0 / std::numbers::pi);

    REQUIRE_THAT(mean, Catch::Matchers::WithinAbs(expected_mean, 1e-2));
    REQUIRE_THAT(var, Catch::Matchers::WithinAbs(expected_var, 1e-2));
}

TEST_CASE(
    "HalfNormalSampler reset clears cached distribution state",
    "[stats][half_normal_sampler]")
{
    HalfNormalSampler<double> sampler{1.0};
    const HalfNormalSampler<double>::Params post{.mean = 0.0, .var = 1.0};

    std::mt19937_64 rng{2};
    (void)sampler.draw(post, +1, rng);
    sampler.reset();
    rng.seed(2);
    const double after_reset = sampler.draw(post, +1, rng);

    HalfNormalSampler<double> fresh{1.0};
    std::mt19937_64 fresh_rng{2};
    const double fresh_first = fresh.draw(post, +1, fresh_rng);

    REQUIRE(after_reset == fresh_first);
}

TEST_CASE(
    "HalfNormalSampler draw uses posterior sign",
    "[stats][half_normal_sampler]")
{
    HalfNormalSampler<double> hn{1.0};
    const NormalSampler<double>::Kernel kernel{
        .quadratic = 1.0,
        .linear = 1.0,
        .scale = 1.0,
    };
    const auto post = hn.posterior_with_logL(kernel, -1);

    std::mt19937_64 rng{SEED};
    for (int i = 0; i < 1000; ++i)
    {
        REQUIRE(hn.draw(post, rng) < 0.0);
    }
}

TEST_CASE(
    "HalfNormalSampler draw truncated normal Devroye path sign=+1 mean=-3",
    "[stats][half_normal_sampler]")
{
    HalfNormalSampler<double> hn{1.0};
    const HalfNormalSampler<double>::Params post{.mean = -3.0, .var = 1.0};

    std::mt19937_64 rng{SEED};
    double sum = 0.0;
    for (int i = 0; i < DRAW_COUNT; ++i)
    {
        sum += hn.draw(post, +1, rng);
    }
    const double sample_mean = sum / DRAW_COUNT;

    constexpr double MU = -3.0;
    constexpr double SIGMA = 1.0;
    constexpr double ALPHA = 3.0;  // (0 - MU) / SIGMA
    const double phi_alpha = std::exp((-0.5 * ALPHA) * ALPHA)
                             / std::sqrt(2.0 * std::numbers::pi);
    const double Phi_alpha = 0.5 * std::erfc(-ALPHA / std::numbers::sqrt2);
    const double expected_mean = MU + (SIGMA * phi_alpha / (1.0 - Phi_alpha));

    REQUIRE_THAT(sample_mean, Catch::Matchers::WithinAbs(expected_mean, 1e-2));
}

}  // namespace
}  // namespace gelex
