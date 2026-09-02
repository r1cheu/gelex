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

#ifndef GELEX_TEST_COMPACT_GENOTYPE_FIXTURE_H_
#define GELEX_TEST_COMPACT_GENOTYPE_FIXTURE_H_

#include <Eigen/Core>
#include <utility>
#include <vector>

#include "gelex/bayes/design.h"
#include "gelex/bayes/model.h"
#include "gelex/data/bed.h"
#include "gelex/data/genotype_method.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"

#include "bed_fixture.h"

namespace gelex::test
{

inline auto make_bed(Eigen::MatrixXd genotypes) -> Bed
{
    BedFixture fixture;
    const auto [prefix, raw]
        = fixture.create_deterministic_bed_files(genotypes);
    static_cast<void>(raw);
    return open_bed(prefix.string());
}

inline auto make_genetic_design(
    Eigen::MatrixXd genotypes,
    GeneticModeSet modes = GeneticModeSet{GeneticMode::A},
    GenotypeMethod method = GenotypeMethod::Center) -> bayes::GeneticDesign
{
    return bayes::GeneticDesign{make_bed(std::move(genotypes)), modes, method};
}

inline auto make_genetic_design_without_modes(Eigen::MatrixXd genotypes)
    -> bayes::GeneticDesign
{
    return bayes::GeneticDesign{make_bed(std::move(genotypes))};
}

inline auto make_compact_model(
    Eigen::MatrixXd genotypes,
    Eigen::VectorXd phenotype,
    GeneticModeSet modes = GeneticModeSet{GeneticMode::A},
    GenotypeMethod method = GenotypeMethod::Center) -> BayesModel
{
    auto genetic = make_genetic_design(genotypes, modes, method);
    return BayesModel{
        std::move(phenotype),
        FixedDesign::make(genotypes.rows()),
        {},
        std::move(genetic)};
}

inline auto make_compact_model(
    Eigen::MatrixXd genotypes,
    GeneticModeSet modes = GeneticModeSet{GeneticMode::A},
    GenotypeMethod method = GenotypeMethod::Center) -> BayesModel
{
    return make_compact_model(
        std::move(genotypes),
        Eigen::VectorXd{{1.0, -0.5, 0.25}},
        modes,
        method);
}

}  // namespace gelex::test

#endif  // GELEX_TEST_COMPACT_GENOTYPE_FIXTURE_H_
