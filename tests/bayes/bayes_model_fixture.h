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

#ifndef GELEX_TEST_BAYES_BAYES_MODEL_FIXTURE_H_
#define GELEX_TEST_BAYES_BAYES_MODEL_FIXTURE_H_

#include <Eigen/Core>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/design_data.h"
#include "gelex/bayes/model.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"

#include "compact_genotype_fixture.h"

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
    random.emplace_back(
        "batch",
        std::vector<std::string>{"batch_1", "batch_2"},
        Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}});
    return BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        FixedDesign::make(3),
        std::move(random),
        std::move(genetic)};
}

}  // namespace gelex::test

#endif  // GELEX_TEST_BAYES_BAYES_MODEL_FIXTURE_H_
