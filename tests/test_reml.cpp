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
#include <catch2/matchers/catch_matchers.hpp>
#include <cmath>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/algo/reml/estimator.h"
#include "gelex/algo/reml/optimizer_state.h"
#include "gelex/algo/reml/variance_calculator.h"  // IWYU pragma: keep
#include "gelex/freq/design.h"
#include "gelex/freq/model.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_effect_type.h"

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
    double sigma_g{};
};

auto make_closed_form_problem() -> RemlProblem
{
    Eigen::VectorXd y(3);
    y << 4.0, 1.0, -1.0;

    Eigen::MatrixXd K(3, 3);
    // clang-format off
    K <<  8.0 / 3.0,  2.0 / 3.0, -1.0 / 3.0,
          2.0 / 3.0,  8.0 / 3.0, -1.0 / 3.0,
         -1.0 / 3.0, -1.0 / 3.0, 11.0 / 3.0;
    // clang-format on

    return {
        .y = y,
        .X = FixedDesign::build(y.size()).X,
        .K = K,
        .sigma_e = 5.0 / 6.0,
        .sigma_g = 11.0 / 6.0};
}

auto build_model(const RemlProblem& problem) -> std::pair<FreqModel, FreqState>
{
    auto fixed = FixedDesign::build(problem.y.size());
    fixed.X = problem.X;

    std::vector<freq::GeneticDesign> genetics;
    genetics.push_back({.type = GeneticMode::A, .K = problem.K});

    FreqModel model(problem.y, std::move(fixed), std::move(genetics));
    FreqState state(model);
    state.residual().variance = problem.sigma_e;
    state.genetic()[0].variance = problem.sigma_g;
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
          + problem.sigma_g * problem.K;
    const Eigen::MatrixXd Vinv_ref = V_ref.inverse();
    const Eigen::MatrixXd XtVinvX_ref
        = problem.X.transpose() * Vinv_ref * problem.X;
    const Eigen::MatrixXd XtVinvX_inv_ref = XtVinvX_ref.inverse();
    const Eigen::MatrixXd P_ref = Vinv_ref
                                  - Vinv_ref * problem.X * XtVinvX_inv_ref
                                        * problem.X.transpose() * Vinv_ref;
    const Eigen::VectorXd Py_ref = P_ref * problem.y;

    reml::OptimizerState opt(model);
    reml::compute_v(model, state, opt.V);
    opt.logdet_v = reml::v_inv_logdet(opt.V);
    reml::compute_proj(model, opt);

    SECTION("in-place materialized P matches closed form")
    {
        Eigen::MatrixXd P_mat = opt.V;
        P_mat.noalias() -= opt.ViX * opt.XtViX_inv * opt.ViX.transpose();
        REQUIRE(P_mat.isApprox(P_ref, 1e-10));
        REQUIRE(opt.Py.isApprox(Py_ref, 1e-10));
    }

    SECTION("compute_loglike: -0.5*(log|V| + log|X'V^{-1}X| + y'Py)")
    {
        const double expected
            = -0.5
              * (std::log(V_ref.determinant())
                 + std::log(XtVinvX_ref.determinant()) + problem.y.dot(Py_ref));
        REQUIRE(
            std::abs(reml::compute_loglike(model, opt) - expected) <= 1e-10);
    }
}

TEST_CASE(
    "Estimator::fit — 3x3 two-component matches REML closed form",
    "[reml][canary]")
{
    // n=3, intercept-only X (p=1), single genetic component. K is built as
    //   K = 3*v0 v0' + 2*u1 u1' + 4*u2 u2',
    //   v0 = [1,1,1]/sqrt(3), u1 = [1,-1,0]/sqrt(2), u2 = [1,1,-2]/sqrt(6)
    // so K_tilde = X_perp' K X_perp = diag(2, 4). With y = [4, 1, -1]' the
    // saturated REML MLE solves  sigma_e + sigma_g*lambda_i = z_i^2  giving
    //   sigma_e = 5/6,  sigma_g = 11/6
    // (cross-checked against HIBLUP --single-trait).
    const auto problem = make_closed_form_problem();
    auto [model, state] = build_model(problem);

    reml::Estimator estimator(/*max_iter=*/500, /*tol=*/1e-12);
    estimator.fit(model, state);

    REQUIRE(estimator.is_converged());
    REQUIRE_THAT(
        state.residual().variance,
        Catch::Matchers::WithinAbs(problem.sigma_e, 1e-8));
    REQUIRE_THAT(
        state.genetic()[0].variance,
        Catch::Matchers::WithinAbs(problem.sigma_g, 1e-8));
}

TEST_CASE("Estimator::fit resets convergence state between runs", "[reml]")
{
    const auto problem = make_closed_form_problem();
    reml::Estimator estimator(/*max_iter=*/500, /*tol=*/1e-12);

    auto [first_model, first_state] = build_model(problem);
    estimator.fit(first_model, first_state);
    const auto first_iter_count = estimator.iter_count();

    auto [second_model, second_state] = build_model(problem);
    estimator.fit(second_model, second_state);

    REQUIRE(estimator.is_converged());
    REQUIRE(first_iter_count > 1);
    REQUIRE(estimator.iter_count() == first_iter_count);
}

}  // namespace gelex
