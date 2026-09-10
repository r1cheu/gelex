// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <string>

#include "gelex/bayes/draws_diagnostics.h"
#include "gelex/bayes/stats/diagnostics.h"
#include "gelex/exception.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/dense_reader.h"
#include "gelex/io/dense_writer.h"

#include "file_fixture.h"

namespace
{

constexpr std::uint64_t n_draws = 64;

auto random_matrix(std::mt19937_64& rng, Eigen::Index rows) -> Eigen::MatrixXd
{
    std::normal_distribution<double> dist;
    return Eigen::MatrixXd::NullaryExpr(
        rows, static_cast<Eigen::Index>(n_draws), [&]() { return dist(rng); });
}

auto write_payload(
    gelex::DenseWriter& writer,
    const std::string& identifier,
    const Eigen::MatrixXd& values) -> void
{
    auto stream = writer.reserve<double>(
        identifier,
        gelex::BinaryShape{static_cast<std::uint64_t>(values.rows()), n_draws});
    stream << std::span<const double>{
        values.data(), static_cast<std::size_t>(values.size())};
}

auto same_diagnostics(
    const gelex::ChainDiagnostics& lhs,
    const gelex::ChainDiagnostics& rhs) -> bool
{
    return lhs.mean == rhs.mean && lhs.sd == rhs.sd && lhs.median == rhs.median
           && lhs.hpdi_lower == rhs.hpdi_lower
           && lhs.hpdi_upper == rhs.hpdi_upper && lhs.ess == rhs.ess
           && lhs.mcse == rhs.mcse && lhs.split_rhat == rhs.split_rhat;
}

}  // namespace

TEST_CASE("Draws diagnostics per model term", "[bayes][draws_diagnostics]")
{
    gelex::test::FileFixture fixture;
    const auto path = (fixture.get_test_dir() / "terms.draws").string();

    std::mt19937_64 rng(11);
    const Eigen::MatrixXd fixed = random_matrix(rng, 3);
    const Eigen::MatrixXd group = random_matrix(rng, 4);
    const Eigen::MatrixXd group_variance = random_matrix(rng, 1).cwiseAbs();
    const Eigen::MatrixXd residual = random_matrix(rng, 1).cwiseAbs();
    const Eigen::MatrixXd wide_variance = random_matrix(rng, 2);

    {
        gelex::DenseWriter writer{path};
        write_payload(writer, "fixed/coefficients", fixed);
        write_payload(writer, "random/Group/coefficients", group);
        write_payload(writer, "random/Group/variance", group_variance);
        write_payload(writer, "random/Wide/coefficients", group);
        write_payload(writer, "random/Wide/variance", wide_variance);
        write_payload(writer, "residual/variance", residual);
        writer.close();
    }
    const gelex::DenseReader reader{path};

    SECTION("fixed coefficients follow payload rows")
    {
        const auto result = gelex::diagnose_fixed(reader, 0.9);
        REQUIRE(result.size() == 3);
        for (Eigen::Index row = 0; row < fixed.rows(); ++row)
        {
            REQUIRE(same_diagnostics(
                result[static_cast<std::size_t>(row)],
                gelex::diagnose_chain(fixed.row(row), 0.9)));
        }
    }

    SECTION("random effect returns levels and variance")
    {
        const auto result = gelex::diagnose_random(reader, "Group");
        REQUIRE(result.coefficients.size() == 4);
        REQUIRE(same_diagnostics(
            result.coefficients[2], gelex::diagnose_chain(group.row(2))));
        REQUIRE(same_diagnostics(
            result.variance, gelex::diagnose_chain(group_variance.row(0))));
    }

    SECTION("residual variance")
    {
        REQUIRE(same_diagnostics(
            gelex::diagnose_residual(reader),
            gelex::diagnose_chain(residual.row(0))));
    }

    SECTION("probability is forwarded to the hpdi")
    {
        const auto full = gelex::diagnose_residual(reader, 1.0);
        REQUIRE(full.hpdi_lower == residual.minCoeff());
        REQUIRE(full.hpdi_upper == residual.maxCoeff());
    }

    SECTION("missing payload throws")
    {
        REQUIRE_THROWS_AS(
            gelex::diagnose_random(reader, "Missing"), gelex::GelexException);
    }

    SECTION("variance payload with several rows throws")
    {
        REQUIRE_THROWS_AS(
            gelex::diagnose_random(reader, "Wide"), gelex::GelexException);
    }
}
