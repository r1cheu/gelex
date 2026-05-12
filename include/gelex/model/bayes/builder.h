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

#include <span>
#include <vector>

#include "gelex/model/bayes/legacy_method.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

struct GeneticStats
{
    GeneticMode mode{};
    double marker_variance_sum{};
    double heritability_init{};
};

auto build_bayes_method(
    const LegacyBayesConfig&,
    std::span<const GeneticStats>,
    double phenotype_variance) -> LegacyBayesMethod;

}  // namespace gelex::bayes

namespace gelex
{

class PhenoPipe;
class GenoPipe;
class BayesModel;

auto build_bayes_model(PhenoPipe&& pheno, GenoPipe&& geno) -> BayesModel;

auto compute_genetic_stats(const BayesModel& model)
    -> std::vector<bayes::GeneticStats>;

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_BUILDER_H_
