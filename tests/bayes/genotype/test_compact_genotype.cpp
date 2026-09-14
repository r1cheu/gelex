// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/genotype/compact_genotype.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/snp_lut.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

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

    std::vector<std::size_t> completed_markers;
    const auto genotype = gelex::bayes::make_compact_genotype(
        bed,
        [&](std::size_t current) { completed_markers.push_back(current); });

    REQUIRE(genotype.rows() == 3);
    REQUIRE(genotype.cols() == 3);
    REQUIRE(genotype.a1_frequency().isApprox(
        Eigen::VectorXd{{1.0 / 6.0, 2.0 / 3.0, 0.0}}));
    REQUIRE(completed_markers == std::vector<std::size_t>{1, 2, 3});

    // Gathered order is d, b, a; dosage 0/1/2 decodes to codes 3/2/0.
    REQUIRE(
        std::vector<std::uint8_t>{
            genotype.col(0).begin(), genotype.col(0).end()}
        == std::vector<std::uint8_t>{3, 2, 3});
    REQUIRE(genotype.locus_stats()[0].nA2A2 == 2);
    REQUIRE(genotype.locus_stats()[2].n_missing == 3);
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
        auto genetic = gelex::bayes::make_genetic_design(
            gelex::open_bed(prefix.string()),
            GeneticMode::A | GeneticMode::D,
            method);
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

    auto joint = gelex::test::make_genetic_design(
        Eigen::MatrixXd{{0.0}, {1.0}, {2.0}}, GeneticMode::A | GeneticMode::D);

    REQUIRE_THROWS_AS(single.projection(GeneticMode::D), gelex::GelexException);
    REQUIRE_NOTHROW(joint.projection(GeneticMode::A));
    REQUIRE_NOTHROW(joint.projection(GeneticMode::D));
}

TEST_CASE(
    "GeneticDesign validates its representation on construction",
    "[bayes][compact]")
{
    using gelex::bayes::GeneticDesign;
    const auto spec = gelex::encoding_spec_from_method(
        GeneticMode::A, GenotypeMethod::Center);
    auto bed = gelex::test::make_bed(Eigen::MatrixXd{{0.0}, {1.0}, {2.0}});
    auto genotype = std::make_unique<gelex::bayes::CompactGenotype>(
        gelex::bayes::make_compact_genotype(bed));
    const auto make_projections = [&](const gelex::bayes::CompactGenotype& g)
    {
        GeneticDesign::projection_array_type projections;
        projections.at(std::to_underlying(GeneticMode::A))
            = gelex::bayes::make_genetic_projection(g, spec);
        return projections;
    };
    const auto make_metadata = []
    {
        return std::move(
                   gelex::test::make_bed(Eigen::MatrixXd{{0.0}, {1.0}, {2.0}}))
            .bim();
    };

    REQUIRE_THROWS_AS(
        GeneticDesign(nullptr, make_projections(*genotype), make_metadata()),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        GeneticDesign(
            std::make_unique<gelex::bayes::CompactGenotype>(
                gelex::bayes::make_compact_genotype(bed)),
            GeneticDesign::projection_array_type{},
            make_metadata()),
        gelex::GelexException);

    const auto stranger = gelex::bayes::make_compact_genotype(bed);
    REQUIRE_THROWS_AS(
        GeneticDesign(
            std::make_unique<gelex::bayes::CompactGenotype>(
                gelex::bayes::make_compact_genotype(bed)),
            make_projections(stranger),
            make_metadata()),
        gelex::GelexException);

    REQUIRE_THROWS_AS(
        GeneticDesign(
            std::make_unique<gelex::bayes::CompactGenotype>(
                gelex::bayes::make_compact_genotype(bed)),
            make_projections(*genotype),
            std::move(
                gelex::test::make_bed(
                    Eigen::MatrixXd{{0.0, 1.0}, {1.0, 2.0}, {2.0, 0.0}}))
                .bim()),
        gelex::GelexException);

    auto projections = make_projections(*genotype);
    const GeneticDesign design{
        std::move(genotype), std::move(projections), make_metadata()};
    REQUIRE(design.rows() == 3);
    REQUIRE(design.cols() == 1);
    REQUIRE(design.contains(GeneticMode::A));
    REQUIRE_FALSE(design.contains(GeneticMode::D));
}

TEST_CASE(
    "GeneticDesign is rebuilt from saved lookup tables",
    "[bayes][compact]")
{
    const auto bed = gelex::test::make_bed(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 2.0}, {2.0, 0.0}});
    const auto trained = gelex::bayes::make_genetic_design(
        bed, GeneticMode::A | GeneticMode::D, GenotypeMethod::NOIACenter);
    gelex::ModeMap<gelex::SnpLutMatrix> luts;
    for (const auto mode : trained.each_mode())
    {
        luts.emplace(mode, trained.projection(mode).snp_luts());
    }

    const auto rebuilt = gelex::bayes::make_genetic_design(bed, luts);
    REQUIRE(rebuilt.rows() == 3);
    REQUIRE(rebuilt.marker_metadata().rows() == 2);
    for (const auto mode : gelex::all_genetic_modes)
    {
        REQUIRE(rebuilt.projection(mode).snp_luts().isApprox(
            trained.projection(mode).snp_luts()));
        REQUIRE(rebuilt.projection(mode).xtx_diag().isApprox(
            trained.projection(mode).xtx_diag()));
    }

    REQUIRE_THROWS_AS(
        gelex::bayes::make_genetic_design(
            bed, gelex::ModeMap<gelex::SnpLutMatrix>{}),
        gelex::GelexException);
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
    const auto design = gelex::bayes::make_genetic_design(
        gelex::open_bed(prefix.string()),
        gelex::GeneticModeSet{GeneticMode::A},
        GenotypeMethod::Center);

    const auto& metadata = design.marker_metadata();
    REQUIRE(
        std::vector<std::string>{
            metadata.index().keys().begin(), metadata.index().keys().end()}
        == std::vector<std::string>{"marker_1", "marker_2"});
    REQUIRE(
        std::vector<std::string>{
            metadata["CHR"].as<std::string>().begin(),
            metadata["CHR"].as<std::string>().end()}
        == std::vector<std::string>{"3", "7"});
    REQUIRE(
        std::vector<std::int32_t>{
            metadata["BP"].as<std::int32_t>().begin(),
            metadata["BP"].as<std::int32_t>().end()}
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
    const auto genotype = gelex::bayes::make_compact_genotype(
        gelex::test::make_bed(Eigen::MatrixXd{{0.0}, {1.0}, {2.0}}));
    REQUIRE(genotype.rows() == 3);
    REQUIRE(genotype.cols() == 1);
}

TEST_CASE(
    "CompactGenotype is constructible from its representation",
    "[bayes][compact]")
{
    using raw_matrix_type = gelex::bayes::CompactGenotype::raw_matrix_type;
    const raw_matrix_type raw_codes{
        {std::uint8_t{0}, std::uint8_t{2}}, {std::uint8_t{3}, std::uint8_t{1}}};
    std::vector<gelex::LocusStats> stats(2);
    stats[0].nA1A1 = 1;
    stats[0].nA2A2 = 1;
    stats[1].nA1A2 = 1;
    stats[1].n_missing = 1;

    const gelex::bayes::CompactGenotype genotype{
        raw_codes, stats, Eigen::VectorXd{{0.5, 0.5}}};
    REQUIRE(genotype.rows() == 2);
    REQUIRE(genotype.cols() == 2);
    REQUIRE(genotype.locus_stats()[1].n_missing == 1);
    REQUIRE(
        std::vector<std::uint8_t>{
            genotype.col(0).begin(), genotype.col(0).end()}
        == std::vector<std::uint8_t>{0, 3});

    REQUIRE_THROWS_AS(
        gelex::bayes::CompactGenotype(raw_codes, stats, Eigen::VectorXd{{0.5}}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        gelex::bayes::CompactGenotype(
            raw_codes,
            std::vector<gelex::LocusStats>(1),
            Eigen::VectorXd{{0.5, 0.5}}),
        gelex::GelexException);
}

TEST_CASE(
    "GeneticProjection is constructible from a genotype and a spec",
    "[bayes][compact]")
{
    const auto spec = gelex::encoding_spec_from_method(
        GeneticMode::A, GenotypeMethod::Center);
    const auto genotype = gelex::bayes::make_compact_genotype(
        gelex::test::make_bed(Eigen::MatrixXd{{0.0}, {1.0}, {2.0}}));
    const auto additive = gelex::bayes::make_genetic_projection(genotype, spec);

    REQUIRE(&additive.genotype() == &genotype);
    Eigen::VectorXd expanded = Eigen::VectorXd::Zero(3);
    additive.axpy(0, 1.0, expanded);
    REQUIRE(expanded.isApprox(Eigen::VectorXd{{-1.0, 0.0, 1.0}}));

    const auto other = gelex::bayes::make_compact_genotype(
        gelex::test::make_bed(Eigen::MatrixXd{{0.0}, {1.0}, {2.0}}));
    REQUIRE_THROWS_AS(
        additive.col_covariance(
            gelex::bayes::make_genetic_projection(other, spec)),
        gelex::GelexException);
}

TEST_CASE(
    "GeneticProjection is constructible from lookup tables",
    "[bayes][compact]")
{
    const auto genotype = gelex::bayes::make_compact_genotype(
        gelex::test::make_bed(
            Eigen::MatrixXd{{0.0, 0.0}, {1.0, 0.0}, {2.0, 0.0}}));
    // Dosage 0/1/2 is stored as code 3/2/0; a centred additive LUT is
    // (1, 0, 0, -1), the monomorphic second marker is left out.
    const gelex::SnpLutMatrix luts{
        {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {-1.0, 0.0}};
    const gelex::bayes::GeneticProjection projection{
        genotype, luts, std::vector<Eigen::Index>{0}};

    REQUIRE(projection.snp_luts().isApprox(luts));
    REQUIRE(projection.xtx_diag().isApprox(Eigen::VectorXd{{2.0, 0.0}}));
    REQUIRE(
        projection.col_var().isApprox(Eigen::RowVectorXd{{2.0 / 3.0, 0.0}}));
    REQUIRE(
        std::vector<Eigen::Index>{
            projection.valid_indices().begin(),
            projection.valid_indices().end()}
        == std::vector<Eigen::Index>{0});

    const auto oracle = gelex::bayes::make_genetic_projection(
        genotype,
        gelex::encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::Center));
    REQUIRE(oracle.snp_luts().isApprox(luts));
    REQUIRE(oracle.xtx_diag().isApprox(projection.xtx_diag()));

    const gelex::bayes::GeneticProjection every_marker{genotype, luts};
    REQUIRE(every_marker.valid_indices().size() == 2);
    REQUIRE(every_marker.xtx_diag().isApprox(projection.xtx_diag()));
    REQUIRE(every_marker.col_var().isApprox(projection.col_var()));

    REQUIRE_THROWS_AS(
        gelex::bayes::GeneticProjection(
            genotype,
            gelex::SnpLutMatrix::Zero(4, 1),
            std::vector<Eigen::Index>{}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        gelex::bayes::GeneticProjection(
            genotype, luts, std::vector<Eigen::Index>{2}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        gelex::bayes::GeneticProjection(
            genotype, luts, std::vector<Eigen::Index>{1, 0}),
        gelex::GelexException);
    REQUIRE_THROWS_AS(
        gelex::bayes::GeneticProjection(
            genotype, luts, std::vector<Eigen::Index>{0, 0}),
        gelex::GelexException);
}
