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

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "bed_fixture.h"
#include "gelex/data/genotype_processor.h"
#include "gelex/data/grm.h"

namespace fs = std::filesystem;

using namespace gelex;  // NOLINT
using Catch::Matchers::WithinAbs;
using gelex::test::are_matrices_equal;
using gelex::test::BedFixture;

TEST_CASE("GRM - Construction with valid BED files", "[grm][construction]")
{
    BedFixture fixture;

    SECTION("Happy path - construct from valid BED prefix")
    {
        const Eigen::Index num_samples = 10;
        const Eigen::Index num_snps = 20;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        REQUIRE_NOTHROW(
            [&]()
            {
                GRM grm(bed_prefix);
                REQUIRE(grm.num_snps() == num_snps);
                REQUIRE(
                    grm.sample_ids().size()
                    == static_cast<size_t>(num_samples));
            }());
    }
}
