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
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/progress.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/genotype_method.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "bayes/genotype/compact_genotype.h"
#include "bed_fixture.h"
#include "compact_genotype_fixture.h"

using gelex::GeneticMode;
using gelex::GenotypeMethod;

TEST_CASE(
    "CompactGenotype caches gathered raw codes and frequencies",
    "[bayes][compact]")
{
    STATIC_REQUIRE(
        !std::is_copy_constructible_v<gelex::bayes::CompactGenotype>);
    STATIC_REQUIRE(std::is_move_constructible_v<gelex::bayes::CompactGenotype>);

    gelex::test::BedFixture fixture;
    const double missing = std::numeric_limits<double>::quiet_NaN();
    const Eigen::MatrixXd genotypes{
        {0.0, 0.0, missing},
        {1.0, 2.0, missing},
        {2.0, 1.0, missing},
        {0.0, 2.0, missing}};
    const std::vector<std::string> sample_ids{"a", "b", "c", "d"};
    const auto prefix
        = fixture.create_deterministic_bed_files(genotypes, sample_ids).first;
    auto bed = gelex::open_bed(prefix.string());
    const auto source_keys = bed.sample_index().keys();
    bed.gather(
        gelex::DataFrameIndex<std::string>{std::vector<std::string>{
            source_keys[3], source_keys[1], source_keys[0]}});

    std::size_t progress_events = 0;
    bool done = false;
    const gelex::bayes::CompactGenotype genotype{
        bed,
        [&](const gelex::GenotypeProgressEvent& event)
        {
            if (event.done)
            {
                done = true;
            }
            else
            {
                ++progress_events;
            }
        }};

    REQUIRE(genotype.rows() == 3);
    REQUIRE(genotype.cols() == 3);
    REQUIRE(genotype.size_bytes() == 9);
    REQUIRE(genotype.a1_frequency().isApprox(
        Eigen::VectorXd{{1.0 / 6.0, 2.0 / 3.0, 0.0}}));
    REQUIRE(progress_events == 3);
    REQUIRE(done);
}

TEST_CASE(
    "Compact designs match dense encodings for every genotype method",
    "[bayes][compact]")
{
    gelex::test::BedFixture fixture;
    const double missing = std::numeric_limits<double>::quiet_NaN();
    const Eigen::MatrixXd genotypes{
        {0.0, 0.0, missing},
        {1.0, 2.0, missing},
        {2.0, 1.0, missing},
        {0.0, 2.0, missing}};
    const auto prefix = fixture.create_deterministic_bed_files(genotypes).first;

    for (const auto method :
         {GenotypeMethod::StandardizeHWE,
          GenotypeMethod::CenterHWE,
          GenotypeMethod::Standardize,
          GenotypeMethod::Center,
          GenotypeMethod::OrthStandardizeHWE,
          GenotypeMethod::OrthCenterHWE,
          GenotypeMethod::OrthStandardize,
          GenotypeMethod::OrthCenter,
          GenotypeMethod::NOIAStandardize,
          GenotypeMethod::NOIACenter})
    {
        auto genetic = gelex::bayes::GeneticDesign{
            gelex::open_bed(prefix.string()),
            GeneticMode::A | GeneticMode::D,
            method};
        const auto oracle_bed = gelex::open_bed(prefix.string());
        const gelex::LocusEncoder encoder{oracle_bed};

        for (const auto mode : gelex::all_genetic_modes)
        {
            const auto& projection = genetic.projection(mode);
            Eigen::MatrixXd dense(genetic.rows(), genetic.cols());
            std::vector<Eigen::Index> valid_indices;
            for (Eigen::Index marker = 0; marker < genetic.cols(); ++marker)
            {
                const auto stats = encoder.count(marker);
                const auto encoding = encoder.encoding(
                    marker,
                    stats,
                    gelex::encoding_spec_from_method(mode, method));
                encoder.expand(marker, encoding, dense.col(marker));
                if (encoding.valid)
                {
                    valid_indices.push_back(marker);
                }

                const Eigen::VectorXd probe{{0.5, -1.0, 2.0, 0.25}};
                Eigen::VectorXd expanded = Eigen::VectorXd::Zero(4);
                projection.axpy(marker, 1.0, expanded);
                CHECK(expanded.isApprox(dense.col(marker)));
                CHECK(
                    projection.dot(marker, probe)
                    == Catch::Approx(dense.col(marker).dot(probe)));
            }

            CHECK(projection.xtx_diag().isApprox(
                dense.colwise().squaredNorm().transpose()));
            Eigen::RowVectorXd variance(dense.cols());
            for (Eigen::Index marker = 0; marker < dense.cols(); ++marker)
            {
                const double mean = dense.col(marker).mean();
                variance[marker]
                    = dense.col(marker).array().square().mean() - mean * mean;
            }
            CHECK(projection.col_var().isApprox(variance));
            CHECK(
                std::vector<Eigen::Index>{
                    projection.valid_indices().begin(),
                    projection.valid_indices().end()}
                == valid_indices);
        }

        const auto& additive_projection = genetic.projection(GeneticMode::A);
        const auto& dominance_projection = genetic.projection(GeneticMode::D);
        Eigen::RowVectorXd covariance(genetic.cols());
        for (Eigen::Index marker = 0; marker < genetic.cols(); ++marker)
        {
            Eigen::VectorXd additive_column = Eigen::VectorXd::Zero(4);
            Eigen::VectorXd dominance_column = Eigen::VectorXd::Zero(4);
            additive_projection.axpy(marker, 1.0, additive_column);
            dominance_projection.axpy(marker, 1.0, dominance_column);
            covariance[marker]
                = (additive_column.array() * dominance_column.array()).mean()
                  - (additive_column.mean() * dominance_column.mean());
        }
        CHECK(additive_projection.col_covariance(dominance_projection)
                  .isApprox(covariance));
        const auto common_valid_indices = genetic.common_valid_indices();
        CHECK(
            std::vector<Eigen::Index>{
                common_valid_indices.begin(), common_valid_indices.end()}
            == std::vector<Eigen::Index>{0, 1});
    }
}

TEST_CASE("GeneticDesign exposes explicit projections", "[bayes][compact]")
{
    auto single = gelex::test::make_genetic_design(
        Eigen::MatrixXd{{0.0}, {1.0}, {2.0}});
    const Eigen::VectorXd probe{{1.0, 2.0, 3.0}};
    Eigen::VectorXd expanded = Eigen::VectorXd::Zero(3);
    const auto& additive = single.projection(GeneticMode::A);

    additive.axpy(0, 1.0, expanded);

    REQUIRE(additive.dot(0, probe) == expanded.dot(probe));
    REQUIRE(additive.xtx_diag().size() == 1);
    REQUIRE(additive.col_var().size() == 1);
    REQUIRE(additive.valid_indices().size() == 1);
    REQUIRE(additive.snp_luts().cols() == 1);

    auto empty = gelex::test::make_genetic_design_without_modes(
        Eigen::MatrixXd{{0.0}, {1.0}, {2.0}});
    auto joint = gelex::test::make_genetic_design(
        Eigen::MatrixXd{{0.0}, {1.0}, {2.0}}, GeneticMode::A | GeneticMode::D);

    REQUIRE_THROWS_AS(empty.projection(GeneticMode::A), gelex::GelexException);
    REQUIRE_THROWS_AS(single.projection(GeneticMode::D), gelex::GelexException);
    REQUIRE_NOTHROW(joint.projection(GeneticMode::A));
    REQUIRE_NOTHROW(joint.projection(GeneticMode::D));
}

TEST_CASE("GeneticDesign retains marker metadata", "[bayes][compact]")
{
    gelex::test::BedFixture fixture;
    const auto prefix
        = fixture
              .create_deterministic_bed_files(
                  Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}},
                  {},
                  {"marker_1", "marker_2"},
                  {"3", "7"},
                  {{'A', 'G'}, {'C', 'T'}})
              .first;
    const gelex::bayes::GeneticDesign design{
        gelex::open_bed(prefix.string()),
        gelex::GeneticModeSet{GeneticMode::A},
        GenotypeMethod::Center};

    const auto& metadata = design.marker_metadata();
    REQUIRE(
        std::vector<std::string>{
            metadata.index().keys().begin(), metadata.index().keys().end()}
        == std::vector<std::string>{"marker_1", "marker_2"});
    REQUIRE(
        std::vector<std::string>{
            metadata["chrom"].as<std::string>().begin(),
            metadata["chrom"].as<std::string>().end()}
        == std::vector<std::string>{"3", "7"});
    REQUIRE(
        std::vector<std::int32_t>{
            metadata["pos"].as<std::int32_t>().begin(),
            metadata["pos"].as<std::int32_t>().end()}
        == std::vector<std::int32_t>{1, 2});
    REQUIRE(
        std::vector<std::string>{
            metadata["A1"].as<std::string>().begin(),
            metadata["A1"].as<std::string>().end()}
        == std::vector<std::string>{"A", "C"});
    REQUIRE(
        std::vector<std::string>{
            metadata["A2"].as<std::string>().begin(),
            metadata["A2"].as<std::string>().end()}
        == std::vector<std::string>{"G", "T"});
}

TEST_CASE("CompactGenotype supports a single marker BED", "[bayes][compact]")
{
    const gelex::bayes::CompactGenotype genotype{
        gelex::test::make_bed(Eigen::MatrixXd{{0.0}, {1.0}, {2.0}})};
    REQUIRE(genotype.rows() == 3);
    REQUIRE(genotype.cols() == 1);
    REQUIRE(genotype.size_bytes() == 3);
}
