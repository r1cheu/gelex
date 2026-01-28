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

#include "../src/types/bayes_effects.h"
#include "../src/utils/math_utils.h"
#include "gelex/data/data_pipe.h"
#include "gelex/data/genotype_mmap.h"

namespace gelex
{
using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::VectorXd;

BayesModel::BayesModel(DataPipe& data_pipe)
    : phenotype_(std::move(data_pipe).take_phenotype())
{
    num_individuals_ = phenotype_.rows();         // NOLINT
    phenotype_var_ = detail::var(phenotype_)(0);  // NOLINT

    auto fixed_effect_names = data_pipe.fixed_effect_names();
    auto fixed_effects = std::move(data_pipe).take_fixed_effects();
    add_fixed_effect(std::move(fixed_effect_names), std::move(fixed_effects).X);

    std::visit(
        [&](auto&& arg) { add_additive(std::forward<decltype(arg)>(arg)); },
        std::move(data_pipe).take_additive_matrix());

    if (data_pipe.has_dominance_matrix())
    {
        std::visit(
            [&](auto&& arg)
            { add_dominance(std::forward<decltype(arg)>(arg)); },
            std::move(data_pipe).take_dominance_matrix());
    }
}

void BayesModel::add_fixed_effect(
    std::vector<std::string>&& levels,
    MatrixXd&& design_matrix)
{
    fixed_.emplace(std::move(levels), std::move(design_matrix));
}

void BayesModel::add_random_effect(
    std::vector<std::string>&& levels,
    MatrixXd&& design_matrix)
{
    random_.emplace_back(std::move(levels), std::move(design_matrix));
}

void BayesModel::add_additive(GenotypeMap&& matrix)
{
    additive_.emplace(std::move(matrix));
}

void BayesModel::add_additive(GenotypeMatrix&& matrix)
{
    additive_.emplace(std::move(matrix));
}

void BayesModel::add_dominance(GenotypeMap&& matrix)
{
    dominant_.emplace(std::move(matrix));
}

void BayesModel::add_dominance(GenotypeMatrix&& matrix)
{
    dominant_.emplace(std::move(matrix));
}

BayesState::BayesState(const BayesModel& model)
{
    if (const auto* effect = model.fixed(); effect)
    {
        fixed_.emplace(*effect);
    }

    if (const auto* effect = model.additive(); effect)
    {
        additive_.emplace(*effect);
    }

    if (const auto* effect = model.dominant(); effect)
    {
        dominant_.emplace(*effect);
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

    if (additive_)
    {
        sum_var += additive_->variance;
    }
    if (dominant_)
    {
        sum_var += dominant_->variance;
    }

    sum_var += residual_.variance;

    if (additive_)
    {
        additive_->heritability = additive_->variance / sum_var;
    }
    if (dominant_)
    {
        dominant_->heritability = dominant_->variance / sum_var;
    }
}

}  // namespace gelex
