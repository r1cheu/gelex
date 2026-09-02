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

#include "gelex/freq/gwas/loco_grm.h"

TEST_CASE("LocoGrmBuilder - Basic Calculation", "[freq][gwas][loco]")
{
    const Eigen::Index n = 3;

    Eigen::MatrixXd x_w = Eigen::MatrixXd::Random(n, 10);
    Eigen::MatrixXd x_i = Eigen::MatrixXd::Random(n, 3);

    Eigen::MatrixXd raw_g_w = x_w * x_w.transpose();
    Eigen::MatrixXd raw_g_i = x_i * x_i.transpose();

    const double k_w = raw_g_w.trace() / static_cast<double>(n);
    const double k_i = raw_g_i.trace() / static_cast<double>(n);

    gelex::LocoGrmBuilder builder(raw_g_w);
    Eigen::MatrixXd loco_grm;
    builder.build_into(raw_g_i, loco_grm);

    const Eigen::MatrixXd expected = (raw_g_w - raw_g_i) / (k_w - k_i);
    REQUIRE(loco_grm.isApprox(expected));
}

TEST_CASE("LocoGrmBuilder - Rejects incompatible shapes", "[freq][gwas][loco]")
{
    const Eigen::MatrixXd whole_grm = Eigen::MatrixXd::Identity(3, 3);
    const Eigen::MatrixXd chromosome_grm = Eigen::MatrixXd::Identity(2, 2);
    gelex::LocoGrmBuilder builder(whole_grm);
    Eigen::MatrixXd target;

    REQUIRE_THROWS(builder.build_into(chromosome_grm, target));
}
