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
#include <catch2/catch_test_macros.hpp>

#include "gelex/bayes/design.h"
#include "gelex/bayes/state.h"
#include "gelex/data/genotype_method.h"
#include "gelex/types/genetic_mode.h"

#include "algo/mcmc/invariant.h"
#include "compact_genotype_fixture.h"

TEST_CASE(
    "Residual adjustment guard maintains the adjusted response",
    "[mcmc][invariant]")
{
    SECTION("dense column")
    {
        const Eigen::VectorXd column{{1.0, -2.0, 0.5}};
        gelex::bayes::ResidualState residual{
            .y_adj = Eigen::VectorXd{{0.5, 1.0, -1.5}}, .variance = 1.0};
        const Eigen::VectorXd initial = residual.y_adj;
        double coefficient = 0.75;

        {
            gelex::DenseResidualAdjustmentGuard guard{
                column, coefficient, residual};
            coefficient = -0.25;
        }

        REQUIRE(residual.y_adj.isApprox(initial + column));
    }

    SECTION("compact genetic projection")
    {
        auto genetic = gelex::test::make_genetic_design(
            Eigen::MatrixXd{{0.0}, {1.0}, {2.0}},
            gelex::GeneticModeSet{gelex::GeneticMode::A},
            gelex::GenotypeMethod::Center);
        Eigen::VectorXd column = Eigen::VectorXd::Zero(genetic.rows());
        genetic.axpy(0, 1.0, column);

        gelex::bayes::ResidualState residual{
            .y_adj = Eigen::VectorXd{{0.5, 1.0, -1.5}}, .variance = 1.0};
        const Eigen::VectorXd initial = residual.y_adj;
        gelex::bayes::GeneticState state{1, genetic.rows()};
        state.coeffs(0) = 0.75;

        {
            gelex::GeneticResidualAdjustmentGuard guard{
                genetic, gelex::GeneticMode::A, 0, state, residual};
            state.coeffs(0) = -0.25;
        }

        REQUIRE(residual.y_adj.isApprox(initial + column));
    }

    SECTION("unchanged coefficient")
    {
        const Eigen::VectorXd column{{1.0, -2.0, 0.5}};
        gelex::bayes::ResidualState residual{
            .y_adj = Eigen::VectorXd{{0.5, 1.0, -1.5}}, .variance = 1.0};
        const Eigen::VectorXd initial = residual.y_adj;
        double coefficient = 0.75;

        {
            gelex::DenseResidualAdjustmentGuard guard{
                column, coefficient, residual};
        }

        REQUIRE(residual.y_adj.isApprox(initial));
    }
}
