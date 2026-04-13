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

#include "gelex/model/bayes/model.h"

#include <string>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <Eigen/Core>

#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/effects.h"

namespace gelex
{
using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::VectorXd;

BayesModel::BayesModel(
    Eigen::VectorXd phenotype,
    FixedEffect fixed_effects,
    std::vector<bayes::GeneticEffect> genetics)
    : phenotype_(std::move(phenotype)), genetics_(std::move(genetics))
{
    num_individuals_ = phenotype_.rows();
    phenotype_var_ = detail::var(phenotype_)(0);
    add_fixed_effect(std::move(fixed_effects));
}

void BayesModel::add_fixed_effect(FixedEffect&& effect)
{
    fixed_ = std::move(effect);
}

void BayesModel::add_random_effect(
    std::string name,
    std::vector<std::string>&& levels,
    MatrixXd&& X)
{
    random_.emplace_back(std::move(name), std::move(levels), std::move(X));
}

}  // namespace gelex
