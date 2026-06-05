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
#include <fstream>
#include <optional>
#include <sstream>
#include <string>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/algo/gwas/assoc_tester.h"
#include "gelex/algo/gwas/assoc_type.h"
#include "gelex/algo/gwas/joint_tester.h"
#include "gelex/algo/reml/result.h"
#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/reader.h"
#include "gelex/data/genotype/process_method.h"
#include "gelex/io/gwas/writer.h"

using gelex::AssocType;
using gelex::GenotypeProcessMethod;
using gelex::JointTester;
using gelex::RemlResult;
using gelex::TestResult;
using gelex::TestResults;
using gelex::dataframe::ColumnType;
using gelex::dataframe::read_dataframe;
using gelex::dataframe::ReadOptions;
using gelex::gwas::GwasWriter;
using gelex::test::FileFixture;

TEST_CASE("JointTester reports df=2 additive-dominance Wald p", "[gwas]")
{
    JointTester tester(GenotypeProcessMethod::Center());
    tester.resize(4, 1);

    auto raw = tester.genotype_buffer();
    raw = Eigen::MatrixXd{{0.0}, {1.0}, {1.0}, {2.0}};

    RemlResult reml;
    reml.P = Eigen::MatrixXd::Identity(4, 4);
    reml.Py = Eigen::VectorXd{{1.0, 4.0, 2.0, 3.0}};
    reml.Vp = 10.0;

    auto results = tester.run(reml);

    const Eigen::VectorXd z_a{{-1.0, 0.0, 0.0, 1.0}};
    const Eigen::VectorXd z_d{{-0.5, 0.5, 0.5, -0.5}};
    const Eigen::Matrix2d XtPX{
        {z_a.dot(z_a), z_a.dot(z_d)}, {z_a.dot(z_d), z_d.dot(z_d)}};
    const Eigen::Vector2d XtPy{{z_a.dot(reml.Py), z_d.dot(reml.Py)}};
    const double det = (XtPX(0, 0) * XtPX(1, 1)) - (XtPX(0, 1) * XtPX(1, 0));
    const Eigen::Matrix2d inv{
        {XtPX(1, 1) / det, -XtPX(0, 1) / det},
        {-XtPX(1, 0) / det, XtPX(0, 0) / det}};
    const Eigen::Vector2d beta = inv * XtPy;
    const double stat_ad = beta.transpose() * XtPX * beta;
    const double expected_p = std::exp(-0.5 * stat_ad);

    REQUIRE(results.joint_p.has_value());
    REQUIRE(std::abs((*results.joint_p)[0] - expected_p) < 1e-12);
}

TEST_CASE("GwasWriter writes joint P_AD between P_D and PVE_D", "[gwas]")
{
    FileFixture files;
    const auto bim_path = files.create_text_file(
        "SNP\tchrom\tpos\tA1\tA2\n"
        "rs1\t1\t123\tA\tG\n",
        ".tsv");

    ReadOptions options;
    options.index_cols = {0};
    constexpr std::array schema{
        ColumnType::String,
        ColumnType::Int,
        ColumnType::String,
        ColumnType::String};
    auto bim = read_dataframe<std::string>(bim_path, options, schema);

    const auto out_prefix = files.generate_random_file_path().string();

    std::array<double, 1> freq{0.25};
    std::array<double, 1> beta_a{1.1};
    std::array<double, 1> se_a{0.2};
    std::array<double, 1> p_a{0.03};
    std::array<double, 1> pve_a{0.04};
    std::array<double, 1> beta_d{5.5};
    std::array<double, 1> se_d{0.6};
    std::array<double, 1> p_d{0.07};
    std::array<double, 1> p_ad{0.08};
    std::array<double, 1> pve_d{0.09};
    std::array<double, 1> pve{0.10};

    TestResults results{
        .freq = std::span<const double>{freq},
        .additive = TestResult{
            .beta = std::span<const double>{beta_a},
            .se = std::span<const double>{se_a},
            .p = std::span<const double>{p_a},
            .pve = std::span<const double>{pve_a},
        },
        .dominance = TestResult{
            .beta = std::span<const double>{beta_d},
            .se = std::span<const double>{se_d},
            .p = std::span<const double>{p_d},
            .pve = std::span<const double>{pve_d},
        },
        .joint_p = std::optional<std::span<const double>>{
            std::span<const double>{p_ad}},
        .total_pve
        = std::optional<std::span<const double>>{std::span<const double>{pve}},
    };

    {
        GwasWriter writer(out_prefix, bim, AssocType::Joint);
        writer.write(0, results);
    }

    std::ifstream ifs(out_prefix + ".gwas.tsv");
    std::ostringstream oss;
    oss << ifs.rdbuf();

    REQUIRE(
        oss.str()
        == "CHR\tSNP\tBP\tA1\tA2\tA1FREQ\t"
           "BETA_A\tSE_A\tP_A\tPVE_A\t"
           "BETA_D\tSE_D\tP_D\tP_AD\tPVE_D\tPVE\n"
           "1\trs1\t123\tA\tG\t0.25\t"
           "1.1\t0.2\t3.000000e-02\t4.000000e-02\t"
           "5.5\t0.6\t7.000000e-02\t8.000000e-02\t"
           "9.000000e-02\t1.000000e-01\n");
}
