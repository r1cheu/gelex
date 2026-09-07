// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_TEST_BAYES_BAYES_MODEL_FIXTURE_H_
#define GELEX_TEST_BAYES_BAYES_MODEL_FIXTURE_H_

#include <Eigen/Core>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/model.h"
#include "gelex/data/fixed_design.h"
#include "gelex/genetic_mode.h"

#include "compact_genotype_fixture.h"
#include "random_design_fixture.h"

namespace gelex::test
{

// Three individuals, two markers and one named two-level random effect: the
// smallest model whose reports carry a row of every kind.
inline auto make_random_effect_model(
    GeneticModeSet modes = GeneticModeSet{GeneticMode::A}) -> BayesModel
{
    auto genetic = make_genetic_design(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}, modes);
    std::vector<bayes::RandomDesign> random;
    random.push_back(make_random_design(
        "batch",
        std::vector<std::string>{"batch_1", "batch_2"},
        Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}}));
    return BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        FixedDesign::make(3),
        std::move(random),
        std::move(genetic)};
}

}  // namespace gelex::test

#endif  // GELEX_TEST_BAYES_BAYES_MODEL_FIXTURE_H_
