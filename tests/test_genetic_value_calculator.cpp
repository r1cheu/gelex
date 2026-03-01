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

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "bed_fixture.h"
#include "gelex/algo/sim/effect_sampler.h"
#include "gelex/algo/sim/genetic_value_calculator.h"

using namespace gelex;  // NOLINT

namespace
{

auto make_effects(Eigen::Index n_snps) -> CausalEffects
{
    CausalEffects effects;
    effects.resize(n_snps);
    effects.additive.setZero();
    effects.dominance.setZero();
    effects.add_class.setZero();
    effects.dom_class.setZero();
    return effects;
}

auto calculate_values(
    const Eigen::Ref<const Eigen::MatrixXd>& genotypes,
    const CausalEffects& effects,
    bool has_dominance) -> GeneticValues
{
    gelex::test::BedFixture bed_fixture;
    auto [bed_path, genotype_matrix]
        = bed_fixture.create_deterministic_bed_files(genotypes);
    (void)genotype_matrix;

    GeneticValueCalculator calculator(bed_path, has_dominance);
    return calculator.calculate(effects);
}

}  // namespace

TEST_CASE("GeneticValueCalculator - additive calculation", "[genetic_calc]")
{
    SECTION("Single SNP with unit effect")
    {
        Eigen::MatrixXd geno(3, 1);
        geno << 0, 1, 2;

        auto effects = make_effects(1);
        effects.additive(0) = 1.0;

        auto result = calculate_values(geno, effects, false);

        REQUIRE(result.additive.size() == 3);
        REQUIRE(result.dominance.size() == 0);
    }

    SECTION("Multiple SNPs with sparse effects")
    {
        Eigen::MatrixXd geno(4, 5);
        geno << 0, 1, 2, 1, 0, 1, 1, 1, 2, 1, 2, 0, 0, 0, 2, 1, 2, 1, 0, 1;

        auto effects = make_effects(5);
        effects.additive(1) = 1.0;
        effects.additive(3) = -0.5;

        auto result = calculate_values(geno, effects, false);

        REQUIRE(result.additive.size() == 4);

        bool all_zero = true;
        for (Eigen::Index i = 0; i < result.additive.size(); ++i)
        {
            if (std::abs(result.additive(i)) > 1e-10)
            {
                all_zero = false;
                break;
            }
        }
        REQUIRE_FALSE(all_zero);
    }
}

TEST_CASE("GeneticValueCalculator - dominance calculation", "[genetic_calc]")
{
    SECTION("Dominance values computed when enabled")
    {
        Eigen::MatrixXd geno(3, 2);
        geno << 0, 2, 1, 1, 2, 0;

        auto effects = make_effects(2);
        effects.additive << 1.0, 0.5;
        effects.dominance << 0.5, 1.0;

        auto result = calculate_values(geno, effects, true);

        REQUIRE(result.additive.size() == 3);
        REQUIRE(result.dominance.size() == 3);

        bool has_nonzero = false;
        for (Eigen::Index i = 0; i < result.dominance.size(); ++i)
        {
            if (std::abs(result.dominance(i)) > 1e-10)
            {
                has_nonzero = true;
                break;
            }
        }
        REQUIRE(has_nonzero);
    }

    SECTION("Dominance vector is empty when disabled")
    {
        Eigen::MatrixXd geno(3, 2);
        geno << 0, 2, 1, 1, 2, 0;

        auto effects = make_effects(2);
        effects.additive << 1.0, 0.5;
        effects.dominance << 0.5, 1.0;

        auto result = calculate_values(geno, effects, false);

        REQUIRE(result.additive.size() == 3);
        REQUIRE(result.dominance.size() == 0);
    }
}
