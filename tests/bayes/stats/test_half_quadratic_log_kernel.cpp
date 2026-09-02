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

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "gelex/bayes/stats/half_quadratic_log_kernel.h"
#include "gelex/bayes/stats/quadratic_log_kernel.h"
#include "gelex/bayes/stats/truncated_normal_distribution.h"
#include "gelex/infra/normal.h"

namespace gelex
{
namespace
{

using Catch::Approx;

template <typename Left, typename Right>
concept Addable
    = requires(const Left& left, const Right& right) { left + right; };

static_assert(Addable<QuadraticLogKernel, HalfQuadraticLogKernel>);
static_assert(Addable<HalfQuadraticLogKernel, QuadraticLogKernel>);
static_assert(!Addable<HalfQuadraticLogKernel, HalfQuadraticLogKernel>);

TEST_CASE(
    "Complementary half-normal prior is normalized on both half lines",
    "[bayes][stats][half_quadratic_log_kernel]")
{
    constexpr double variance = 2.5;
    const auto evaluation = make_half_normal_prior(variance).evaluate();

    for (const HalfLine support : {HalfLine::Negative, HalfLine::Positive})
    {
        const auto parameters = evaluation.truncated_normal_parameters(support);

        REQUIRE(evaluation.log_integral(support) == Approx(0.0).margin(1e-12));
        REQUIRE(parameters.mean() == 0.0);
        REQUIRE(parameters.stddev() == Approx(std::sqrt(variance)));
        REQUIRE(parameters.support() == support);
    }
}

TEST_CASE(
    "Complementary half-line kernel addition shares posterior parameters",
    "[bayes][stats][half_quadratic_log_kernel]")
{
    constexpr double prior_variance = 3.0;
    const QuadraticLogKernel likelihood{4.0, 2.0, 0.0};
    const auto normal_posterior
        = likelihood + make_normal_prior(prior_variance);
    const auto normal_parameters = normal_posterior.normal_parameters();

    const auto prior = make_half_normal_prior(prior_variance);
    const auto posterior = (likelihood + prior).evaluate();
    const auto reverse_posterior = (prior + likelihood).evaluate();

    const auto positive_parameters
        = posterior.truncated_normal_parameters(HalfLine::Positive);
    const auto negative_parameters
        = posterior.truncated_normal_parameters(HalfLine::Negative);
    REQUIRE(positive_parameters.mean() == Approx(normal_parameters.mean()));
    REQUIRE(positive_parameters.stddev() == Approx(normal_parameters.stddev()));
    REQUIRE(positive_parameters.support() == HalfLine::Positive);
    REQUIRE(negative_parameters.mean() == Approx(normal_parameters.mean()));
    REQUIRE(negative_parameters.stddev() == Approx(normal_parameters.stddev()));
    REQUIRE(negative_parameters.support() == HalfLine::Negative);
    REQUIRE(
        reverse_posterior.log_integral(HalfLine::Negative)
        == Approx(posterior.log_integral(HalfLine::Negative)));
    REQUIRE(
        reverse_posterior.log_integral(HalfLine::Positive)
        == Approx(posterior.log_integral(HalfLine::Positive)));
}

TEST_CASE(
    "Complementary half-line integrals add both Gaussian tail masses",
    "[bayes][stats][half_quadratic_log_kernel]")
{
    constexpr double prior_variance = 3.0;
    const QuadraticLogKernel likelihood{4.0, 2.0, 0.0};
    const auto normal_posterior
        = likelihood + make_normal_prior(prior_variance);
    const auto normal_parameters = normal_posterior.normal_parameters();
    const double standardized_mean
        = normal_parameters.mean() / normal_parameters.stddev();
    const double half_normal_normalization = std::log(2.0);

    const auto evaluation
        = (likelihood + make_half_normal_prior(prior_variance)).evaluate();

    REQUIRE(
        evaluation.log_integral(HalfLine::Positive)
        == Approx(
               normal_posterior.log_integral() + half_normal_normalization
               + log_norm_cdf(standardized_mean))
               .margin(1e-12));
    REQUIRE(
        evaluation.log_integral(HalfLine::Negative)
        == Approx(
               normal_posterior.log_integral() + half_normal_normalization
               + log_norm_cdf(-standardized_mean))
               .margin(1e-12));
}

}  // namespace
}  // namespace gelex
