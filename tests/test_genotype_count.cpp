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
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

#include "gelex/data/detail/bed_source.h"
#include "gelex/data/encode/detail/count.h"
#include "gelex/data/encode/detail/sample_mask.h"
#include "gelex/data/encode/stats.h"

#include "bed_fixture.h"
#include "locus_stats_oracle.h"

namespace
{
using gelex::LocusStats;
using gelex::detail::count_genotypes;
using gelex::detail::SampleMask;

// Straightforward per-sample tabulation used as the oracle for the bit-plane
// kernel. Iterates the kept source samples in target order and mirrors the
// decode LUT mapping: 00->nA1A1, 01->missing, 10->nA1A2, 11->nA2A2.
auto reference_count(
    std::span<const std::uint8_t> bytes,
    std::span<const Eigen::Index> target_to_source) -> LocusStats
{
    LocusStats stats;

    for (const Eigen::Index source : target_to_source)
    {
        const auto byte = bytes[static_cast<std::size_t>(source / 4)];
        const auto shift = static_cast<unsigned>(2 * (source % 4));
        const auto code = (static_cast<unsigned>(byte) >> shift) & 0x03U;

        switch (code)
        {
            case 0:
                ++stats.nA1A1;
                break;
            case 1:
                ++stats.n_missing;
                break;
            case 2:
                ++stats.nA1A2;
                break;
            case 3:
                ++stats.nA2A2;
                break;
            default:
                break;
        }
    }

    return stats;
}

auto expect_equal(const LocusStats& got, const LocusStats& want) -> void
{
    CHECK(got.nA1A1 == want.nA1A1);
    CHECK(got.nA1A2 == want.nA1A2);
    CHECK(got.nA2A2 == want.nA2A2);
    CHECK(got.n_missing == want.n_missing);
}

auto identity_projection(Eigen::Index n) -> std::vector<Eigen::Index>
{
    std::vector<Eigen::Index> p(static_cast<std::size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i)
    {
        p[static_cast<std::size_t>(i)] = i;
    }
    return p;
}
}  // namespace

TEST_CASE(
    "count_genotypes matches per-class tabulation",
    "[data][genotype_count]")
{
    SECTION("one full byte, all classes present")
    {
        // samples 0..3 -> codes 00, 01, 10, 11 (little-endian within byte)
        const std::vector<std::uint8_t> bytes{0b11100100};
        const auto proj = identity_projection(4);
        const SampleMask mask{proj, 4};

        const auto stats = count_genotypes(bytes, mask);
        expect_equal(stats, reference_count(bytes, proj));
        CHECK(stats.nA1A1 == 1);
        CHECK(stats.n_missing == 1);
        CHECK(stats.nA1A2 == 1);
        CHECK(stats.nA2A2 == 1);
    }

    SECTION("sample-count padding: source_size not a multiple of 4")
    {
        // 5 samples -> 2 bytes; the top 6 bits of byte 1 are padding and must
        // be ignored even when non-zero.
        const std::vector<std::uint8_t> bytes{0b11100100, 0b11111111};
        const auto proj = identity_projection(5);
        const SampleMask mask{proj, 5};

        const auto stats = count_genotypes(bytes, mask);
        expect_equal(stats, reference_count(bytes, proj));
        CHECK(stats.nA1A1 + stats.nA1A2 + stats.nA2A2 + stats.n_missing == 5);
    }

    SECTION("sample subset keeps only the projected source rows")
    {
        const std::vector<std::uint8_t> bytes{0b11100100};
        // target order keeps source 0 and source 2 out of 4.
        const std::vector<Eigen::Index> proj{0, 2};
        const SampleMask mask{proj, 4};

        const auto stats = count_genotypes(bytes, mask);
        expect_equal(stats, reference_count(bytes, proj));
        CHECK(mask.n_kept() == 2);
    }
}

TEST_CASE(
    "count_genotypes randomized oracle cross-check",
    "[data][genotype_count]")
{
    std::mt19937 rng{0xC0FFEE};
    std::bernoulli_distribution keep_dist{0.7};

    for (const Eigen::Index source_size :
         {Eigen::Index{1},
          Eigen::Index{3},
          Eigen::Index{4},
          Eigen::Index{31},
          Eigen::Index{32},
          Eigen::Index{33},
          Eigen::Index{100},
          Eigen::Index{257}})
    {
        const auto num_bytes = static_cast<std::size_t>((source_size + 3) / 4);
        std::uniform_int_distribution<int> byte_dist{0, 255};

        for (int trial = 0; trial < 32; ++trial)
        {
            std::vector<std::uint8_t> bytes(num_bytes);
            for (auto& b : bytes)
            {
                b = static_cast<std::uint8_t>(byte_dist(rng));
            }

            std::vector<Eigen::Index> proj;
            for (Eigen::Index s = 0; s < source_size; ++s)
            {
                if (keep_dist(rng))
                {
                    proj.push_back(s);
                }
            }

            const SampleMask mask{proj, source_size};
            const auto got = count_genotypes(bytes, mask);
            const auto want = reference_count(bytes, proj);

            INFO("source_size=" << source_size << " trial=" << trial);
            expect_equal(got, want);
            CHECK(mask.n_kept() == static_cast<Eigen::Index>(proj.size()));
        }
    }
}

// End-to-end anchor: pins the kernel to the fixture genotype truth (the source
// that wrote the .bed) by cross-checking against compute_locus_stats over the
// packed bytes, rather than a second hand-written oracle.
TEST_CASE(
    "count_genotypes agrees with the fixture genotype truth",
    "[data][genotype_count]")
{
    gelex::test::BedFixture fixture;
    const auto [prefix, genotypes] = fixture.create_bed_files(
        /*num_samples=*/37,
        /*num_snps=*/24,
        /*missing_rate=*/0.1,
        /*maf_min=*/0.05,
        /*maf_max=*/0.5,
        /*seed=*/42);

    const auto source = gelex::detail::open_bed_source(prefix.string());

    const auto proj = identity_projection(genotypes.rows());
    const SampleMask mask{proj, genotypes.rows()};

    for (Eigen::Index snp = 0; snp < genotypes.cols(); ++snp)
    {
        const auto want
            = gelex::test::compute_locus_stats<double>(genotypes.col(snp));
        const auto got
            = count_genotypes(source[static_cast<std::size_t>(snp)], mask);

        INFO("snp=" << snp);
        expect_equal(got, want);
    }
}
