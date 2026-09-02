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

#include "gelex/bayes/stats/scaled_inv_chi2_log_kernel.h"

using Catch::Approx;
using gelex::make_normal_variance_likelihood;
using gelex::make_scaled_inv_chi2_prior;

TEST_CASE(
    "ScaledInvChi2LogKernel combines a prior with a variance likelihood",
    "[bayes][stats][scaled_inv_chi2_log_kernel]")
{
    const auto posterior = make_scaled_inv_chi2_prior(4.0, 1.0)
                           + make_normal_variance_likelihood(20, 12.5);
    const auto parameters = posterior.scaled_inv_chi2_parameters();

    REQUIRE(parameters.degrees_of_freedom() == Approx(24.0));
    REQUIRE(parameters.scale() == Approx(16.5 / 24.0));
}

TEST_CASE(
    "An empty variance likelihood leaves prior parameters unchanged",
    "[bayes][stats][scaled_inv_chi2_log_kernel]")
{
    const auto posterior = make_scaled_inv_chi2_prior(3.0, 0.5)
                           + make_normal_variance_likelihood(0, 0.0);
    const auto parameters = posterior.scaled_inv_chi2_parameters();

    REQUIRE(parameters.degrees_of_freedom() == Approx(3.0));
    REQUIRE(parameters.scale() == Approx(0.5));
}
