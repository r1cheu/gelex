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
#include <numeric>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/detail/bed_source.h"
#include "gelex/data/encode/detail/count.h"
#include "gelex/data/encode/detail/encoding.h"
#include "gelex/data/encode/detail/expand.h"
#include "gelex/data/encode/detail/sample_mask.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"

#include "bed_fixture.h"
#include "locus_stats_oracle.h"

namespace
{
constexpr std::array<Eigen::Index, 3> raw_code_by_dosage{3, 2, 0};

// Mirrors GrmBuilder's fused path: tabulate once per variant, then expand each
// packed column straight into an encoded column.
auto fused_encode(
    const gelex::Bed& bed,
    const gelex::detail::BedSource& source,
    gelex::GeneticMode mode,
    gelex::GenotypeMethod method) -> Eigen::MatrixXd
{
    const auto num_samples = bed.num_samples();
    const auto num_snps = bed.num_snps();
    std::vector<Eigen::Index> target_to_source(
        static_cast<std::size_t>(num_samples));
    std::iota(
        target_to_source.begin(), target_to_source.end(), Eigen::Index{0});
    const gelex::detail::SampleMask mask{target_to_source, num_samples};
    const auto spec = gelex::encoding_spec_from_method(mode, method);

    Eigen::MatrixXd z(num_samples, num_snps);
    for (Eigen::Index snp = 0; snp < num_snps; ++snp)
    {
        const auto packed = source[static_cast<std::size_t>(snp)];
        const auto stats = gelex::detail::count_genotypes(packed, mask);
        const auto encoding
            = gelex::detail::make_locus_encoding(snp, stats, spec);
        gelex::detail::expand_encoded_column(
            packed, target_to_source, encoding.lut, z.col(snp));
    }
    return z;
}

// Independent dosage oracle for the packed bit-LUT expansion.
auto dosage_encode(
    const Eigen::MatrixXd& dosage,
    gelex::GeneticMode mode,
    gelex::GenotypeMethod method) -> Eigen::MatrixXd
{
    const auto spec = gelex::encoding_spec_from_method(mode, method);
    Eigen::MatrixXd out(dosage.rows(), dosage.cols());
    for (Eigen::Index snp = 0; snp < dosage.cols(); ++snp)
    {
        const auto stats
            = gelex::test::compute_locus_stats<double>(dosage.col(snp));
        const auto encoding
            = gelex::detail::make_locus_encoding(snp, stats, spec);
        for (Eigen::Index i = 0; i < dosage.rows(); ++i)
        {
            const double d = dosage(i, snp);
            const Eigen::Index raw_code
                = std::isnan(d)
                      ? 1
                      : raw_code_by_dosage[static_cast<std::size_t>(d)];
            out(i, snp) = encoding.lut[raw_code];
        }
    }
    return out;
}
}  // namespace

TEST_CASE(
    "expand_encoded_column matches a dense dosage apply",
    "[data][genotype_expand]")
{
    gelex::test::BedFixture fixture;
    const auto [prefix, genotypes] = fixture.create_bed_files(
        /*num_samples=*/41,
        /*num_snps=*/20,
        /*missing_rate=*/0.1,
        /*maf_min=*/0.05,
        /*maf_max=*/0.5,
        /*seed=*/11);

    const auto bed = gelex::open_bed(prefix.string());
    const auto source = gelex::detail::open_bed_source(prefix.string());
    const auto method = gelex::GenotypeMethod::OrthStandardize;

    for (const auto mode : {gelex::GeneticMode::A, gelex::GeneticMode::D})
    {
        const Eigen::MatrixXd expected = dosage_encode(genotypes, mode, method);
        const Eigen::MatrixXd fused = fused_encode(bed, source, mode, method);

        INFO("mode=" << static_cast<int>(mode));
        CHECK(fused.isApprox(expected));
    }
}
