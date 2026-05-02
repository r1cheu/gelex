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

#ifndef GELEX_MODEL_BAYES_BUILDER_H_
#define GELEX_MODEL_BAYES_BUILDER_H_

#include <optional>
#include <vector>

#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/prior.h"

namespace gelex
{

class PhenoPipe;
class GenoPipe;
class BayesModel;

auto build_bayes_model(PhenoPipe&& pheno, GenoPipe&& geno) -> BayesModel;

struct PriorOverrides
{
    BayesMethodConfig method;
    double phenotype_variance{};
    std::optional<std::vector<double>> pi;
    std::optional<std::vector<double>> dpi;
    std::optional<std::vector<double>> multiplier;
    std::optional<std::vector<double>> dmultiplier;
    double positive_prob{0.5};
};

auto build_bayes_priors(
    const PriorOverrides& overrides,
    const BayesModel& model) -> bayes::Priors;

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_BUILDER_H_
