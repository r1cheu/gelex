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

#include <optional>
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
    bayes::GenotypeStorage additive,
    std::optional<bayes::GenotypeStorage> dominance)
    : phenotype_(std::move(phenotype))
{
    num_individuals_ = phenotype_.rows();
    phenotype_var_ = detail::var(phenotype_)(0);
    add_fixed_effect(std::move(fixed_effects));
    genetics_.emplace_back(GeneticKind::Add, std::move(additive));
    if (dominance)
    {
        genetics_.emplace_back(GeneticKind::Dom, std::move(*dominance));
    }
}

void BayesModel::add_fixed_effect(FixedEffect&& effect)
{
    fixed_ = std::move(effect);
}

void BayesModel::add_random_effect(
    std::vector<std::string>&& levels,
    MatrixXd&& X)
{
    random_.emplace_back(std::move(levels), std::move(X));
}

BayesState::BayesState(const BayesModel& model) : fixed_(model.fixed())
{
    for (const auto& effect : model.genetics())
    {
        genetics_.emplace_back(effect);
    }

    if (const auto& effects = model.random(); !effects.empty())
    {
        random_.reserve(effects.size());
        for (const auto& effect : effects)
        {
            random_.emplace_back(effect);
        }
    }
    residual_.y_adj = model.phenotype().array();
    residual_.variance = model.residual().init_variance;
}

void BayesState::compute_heritability()
{
    double sum_var = 0;

    for (const auto& rand : random_)
    {
        sum_var += rand.variance;
    }

    for (const auto& gen : genetics_)
    {
        sum_var += gen.variance;
    }

    sum_var += residual_.variance;

    for (auto& gen : genetics_)
    {
        gen.heritability = gen.variance / sum_var;
    }
}

}  // namespace gelex
