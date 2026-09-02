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
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <utility>

#include "gelex/data/bed.h"
#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/data/genotype_method.h"
#include "gelex/types/genetic_mode.h"

#include "bed_fixture.h"
#include "locus_stats_oracle.h"

namespace
{
constexpr std::array<Eigen::Index, 3> RAW_CODE_BY_DOSAGE{3, 2, 0};

using gelex::LocusEncoder;

constexpr auto METHOD = gelex::GenotypeMethod::OrthStandardize;

struct Dataset
{
    gelex::test::BedFixture fixture;
    std::filesystem::path prefix;
    Eigen::MatrixXd genotypes;
    gelex::Bed bed;
};

auto make_dataset() -> Dataset
{
    gelex::test::BedFixture fixture;
    auto [prefix, genotypes] = fixture.create_bed_files(
        /*num_samples=*/41,
        /*num_snps=*/20,
        /*missing_rate=*/0.1,
        /*maf_min=*/0.05,
        /*maf_max=*/0.5,
        /*seed=*/11);
    auto bed = gelex::open_bed(prefix.string());
    return Dataset{
        std::move(fixture),
        std::move(prefix),
        std::move(genotypes),
        std::move(bed)};
}

// Independent dosage-indexed oracle for the raw BED lookup.
auto expected_column(
    const gelex::LocusEncoding& encoding,
    const Eigen::Ref<const Eigen::VectorXd>& dosage) -> Eigen::VectorXd
{
    Eigen::VectorXd want(dosage.size());
    for (Eigen::Index i = 0; i < dosage.size(); ++i)
    {
        const double d = dosage[i];
        const Eigen::Index raw_code
            = std::isnan(d) ? 1
                            : RAW_CODE_BY_DOSAGE[static_cast<std::size_t>(d)];
        want[i] = encoding.lut[raw_code];
    }
    return want;
}
}  // namespace

TEST_CASE(
    "LocusEncoder::count matches the dosage tabulation oracle",
    "[data][locus_encoder]")
{
    const Dataset data = make_dataset();
    const LocusEncoder encoder{data.bed};

    for (Eigen::Index snp = 0; snp < data.bed.num_snps(); ++snp)
    {
        const gelex::LocusStats got = encoder.count(snp);
        const gelex::LocusStats want
            = gelex::test::compute_locus_stats<double>(data.genotypes.col(snp));

        INFO("snp=" << snp);
        CHECK(got.nA1A1 == want.nA1A1);
        CHECK(got.nA1A2 == want.nA1A2);
        CHECK(got.nA2A2 == want.nA2A2);
        CHECK(got.n_missing == want.n_missing);
    }
}

TEST_CASE(
    "LocusEncoder expand maps each raw code through the LUT per mode",
    "[data][locus_encoder]")
{
    const Dataset data = make_dataset();
    const LocusEncoder encoder{data.bed};
    const auto n = data.bed.num_samples();

    for (const auto mode : {gelex::GeneticMode::A, gelex::GeneticMode::D})
    {
        const auto spec = gelex::encoding_spec_from_method(mode, METHOD);

        for (Eigen::Index snp = 0; snp < data.bed.num_snps(); ++snp)
        {
            const gelex::LocusStats stats = encoder.count(snp);
            const gelex::LocusEncoding encoding
                = encoder.encoding(snp, stats, spec);
            CHECK(encoding.marker_index == snp);

            Eigen::VectorXd got(n);
            encoder.expand(snp, encoding, got);

            INFO("mode=" << static_cast<int>(mode) << " snp=" << snp);
            CHECK(got.isApprox(
                expected_column(encoding, data.genotypes.col(snp))));
        }
    }
}

TEST_CASE("LocusEncoder shares one count across specs", "[data][locus_encoder]")
{
    const Dataset data = make_dataset();
    const LocusEncoder encoder{data.bed};
    const auto n = data.bed.num_samples();

    const auto spec_a
        = gelex::encoding_spec_from_method(gelex::GeneticMode::A, METHOD);
    const auto spec_d
        = gelex::encoding_spec_from_method(gelex::GeneticMode::D, METHOD);

    for (Eigen::Index snp = 0; snp < data.bed.num_snps(); ++snp)
    {
        const gelex::LocusStats stats = encoder.count(snp);
        const gelex::LocusEncoding encoding_a
            = encoder.encoding(snp, stats, spec_a);
        const gelex::LocusEncoding encoding_d
            = encoder.encoding(snp, stats, spec_d);

        Eigen::VectorXd got_a(n);
        Eigen::VectorXd got_d(n);
        encoder.expand(snp, encoding_a, got_a);
        encoder.expand(snp, encoding_d, got_d);

        INFO("snp=" << snp);
        CHECK(got_a.isApprox(
            expected_column(encoding_a, data.genotypes.col(snp))));
        CHECK(got_d.isApprox(
            expected_column(encoding_d, data.genotypes.col(snp))));
    }
}

// Anchors the prediction allele-flip strategy: swapping LUT rows 0/3 is
// equivalent to standardizing 2 - dosage.
TEST_CASE(
    "LocusEncoder::expand with swapped LUT entries equals a 2 - dosage flip",
    "[data][locus_encoder]")
{
    const Dataset data = make_dataset();
    const LocusEncoder encoder{data.bed};
    const auto n = data.bed.num_samples();
    const auto spec
        = gelex::encoding_spec_from_method(gelex::GeneticMode::A, METHOD);

    for (Eigen::Index snp = 0; snp < data.bed.num_snps(); ++snp)
    {
        const gelex::LocusStats stats = encoder.count(snp);
        const gelex::LocusEncoding encoding
            = encoder.encoding(snp, stats, spec);

        gelex::LocusEncoding flipped = encoding;
        std::swap(flipped.lut(0), flipped.lut(3));

        Eigen::VectorXd got(n);
        encoder.expand(snp, flipped, got);

        Eigen::VectorXd want(n);
        for (Eigen::Index i = 0; i < n; ++i)
        {
            const double d = data.genotypes(i, snp);
            if (std::isnan(d))
            {
                want[i] = encoding.lut[1];
                continue;
            }
            const auto flipped_dosage
                = static_cast<std::size_t>(2 - static_cast<Eigen::Index>(d));
            want[i] = encoding.lut[RAW_CODE_BY_DOSAGE[flipped_dosage]];
        }

        INFO("snp=" << snp);
        CHECK(got.isApprox(want));
    }
}
