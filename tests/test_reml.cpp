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
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

#include "gelex/algo/reml/estimator.h"
#include "gelex/algo/reml/reml_buffer.h"
#include "gelex/algo/reml/variance_calculator.h"  // IWYU pragma: keep
#include "gelex/freq/design.h"
#include "gelex/freq/model.h"
#include "gelex/types/fixed_designs.h"

namespace gelex
{
namespace
{

struct RemlProblem
{
    Eigen::VectorXd y;
    Eigen::MatrixXd X;
    Eigen::MatrixXd K;
    double sigma_e{};
    double sigma_random{};
};

auto make_closed_form_problem() -> RemlProblem
{
    Eigen::VectorXd y{{4.0, 1.0, -1.0}};

    Eigen::MatrixXd K{
        {8.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0},
        {2.0 / 3.0, 8.0 / 3.0, -1.0 / 3.0},
        {-1.0 / 3.0, -1.0 / 3.0, 11.0 / 3.0}};

    return {
        .y = y,
        .X = FixedDesign::make(y.size()).X,
        .K = K,
        .sigma_e = 5.0 / 6.0,
        .sigma_random = 11.0 / 6.0};
}

auto build_model(const RemlProblem& problem) -> std::pair<FreqModel, FreqState>
{
    auto fixed = FixedDesign::make(problem.y.size());
    fixed.X = problem.X;

    std::vector<freq::RandomDesign> random;
    random.push_back({.name = "grm", .K = problem.K});

    FreqModel model(problem.y, std::move(fixed), std::move(random));
    FreqState state(model);
    state.residual().variance = problem.sigma_e;
    state.random()[0].variance = problem.sigma_random;
    return {std::move(model), std::move(state)};
}

}  // namespace

TEST_CASE(
    "REML projection — materialized P and loglike at fixed sigma",
    "[reml][canary]")
{
    const auto problem = make_closed_form_problem();
    auto [model, state] = build_model(problem);

    // Reference values via Eigen's generic dense inverse.
    const Eigen::MatrixXd V_ref
        = problem.sigma_e * Eigen::MatrixXd::Identity(3, 3)
          + problem.sigma_random * problem.K;
    const Eigen::MatrixXd Vinv_ref = V_ref.inverse();
    const Eigen::MatrixXd XtVinvX_ref
        = problem.X.transpose() * Vinv_ref * problem.X;
    const Eigen::MatrixXd XtVinvX_inv_ref = XtVinvX_ref.inverse();
    const Eigen::MatrixXd P_ref = Vinv_ref
                                  - Vinv_ref * problem.X * XtVinvX_inv_ref
                                        * problem.X.transpose() * Vinv_ref;
    const Eigen::VectorXd Py_ref = P_ref * problem.y;

    RemlBuffer buffer(model);
    compute_v(model, state, buffer.V);
    buffer.logdet_v = v_inv_logdet(buffer.V);
    compute_proj(model, buffer);

    SECTION("in-place materialized P matches closed form")
    {
        Eigen::MatrixXd P_mat = buffer.V;
        P_mat.noalias()
            -= buffer.ViX * buffer.XtViX_inv * buffer.ViX.transpose();
        REQUIRE(P_mat.isApprox(P_ref, 1e-10));
        REQUIRE(buffer.Py.isApprox(Py_ref, 1e-10));
    }

    SECTION(
        "compute_loglike: -0.5*((n-p)*log(2pi) + log|V| + log|X'V^{-1}X| + "
        "y'Py)")
    {
        const auto dof = static_cast<double>(
            model.num_individuals() - model.fixed().X.cols());
        const double expected
            = -0.5
              * (dof * std::log(2.0 * std::numbers::pi)
                 + std::log(V_ref.determinant())
                 + std::log(XtVinvX_ref.determinant()) + problem.y.dot(Py_ref));
        REQUIRE(std::abs(compute_loglike(model, buffer) - expected) <= 1e-10);
    }
}

TEST_CASE(
    "Estimator::fit — 3x3 two-component matches REML closed form",
    "[reml][canary]")
{
    // n=3, intercept-only X (p=1), single genetic component. K is built as
    // K = 3*v0 v0' + 2*u1 u1' + 4*u2 u2',
    // v0 = [1,1,1]/sqrt(3), u1 = [1,-1,0]/sqrt(2), u2 = [1,1,-2]/sqrt(6)
    // so K_tilde = X_perp' K X_perp = diag(2, 4). With y = [4, 1, -1]' the
    // saturated REML MLE solves sigma_e + sigma_r*lambda_i = z_i^2 giving
    // sigma_e = 5/6, sigma_r = 11/6
    // (cross-checked against HIBLUP --single-trait).
    const auto problem = make_closed_form_problem();
    auto [model, state] = build_model(problem);

    Estimator estimator(/*max_iter=*/500, /*tol=*/1e-12);
    auto fit = estimator.fit(model, state);

    REQUIRE(fit.summary.converged);
    REQUIRE_THAT(
        state.residual().variance,
        Catch::Matchers::WithinAbs(problem.sigma_e, 1e-8));
    REQUIRE_THAT(
        state.random()[0].variance,
        Catch::Matchers::WithinAbs(problem.sigma_random, 1e-8));
}

TEST_CASE("Estimator::fit is deterministic across reuse", "[reml]")
{
    const auto problem = make_closed_form_problem();
    Estimator estimator(/*max_iter=*/500, /*tol=*/1e-12);

    auto [first_model, first_state] = build_model(problem);
    const auto first = estimator.fit(first_model, first_state);

    auto [second_model, second_state] = build_model(problem);
    const auto second = estimator.fit(second_model, second_state);

    REQUIRE(second.summary.converged);
    REQUIRE(first.summary.iter_count > 1);
    REQUIRE(second.summary.iter_count == first.summary.iter_count);
}

}  // namespace gelex
