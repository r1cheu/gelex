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

#include "gelex/bayes/model.h"

#include <algorithm>
#include <ranges>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/bayes/design.h"
#include "gelex/exception.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{
BayesModel::BayesModel(
    Eigen::VectorXd phenotype,
    FixedDesign fixed_design,
    std::vector<bayes::RandomDesign> random,
    std::vector<bayes::GeneticDesign> genetics)
    : phenotype_(std::move(phenotype)),
      fixed_(std::move(fixed_design)),
      random_(std::move(random)),
      genetics_(std::move(genetics))
{
    num_individuals_ = phenotype_.rows();
    phenotype_var_ = stats::detail::vecvar(
        phenotype_, stats::detail::VarNormType::Population);

    if (fixed_.X.rows() != num_individuals_)
    {
        throw GelexException(
            fmt::format(
                "BayesModel: fixed design rows {} != phenotype rows {}",
                fixed_.X.rows(),
                num_individuals_));
    }

    for (const auto& random : random_)
    {
        if (random.X.rows() != num_individuals_)
        {
            throw GelexException(
                fmt::format(
                    "BayesModel: random design '{}' rows {} != phenotype rows "
                    "{}",
                    random.name,
                    random.X.rows(),
                    num_individuals_));
        }
    }

    std::vector<GeneticMode> seen_modes;
    seen_modes.reserve(genetics_.size());
    for (const auto& genetic : genetics_)
    {
        if (genetic.X.rows() != num_individuals_)
        {
            throw GelexException(
                fmt::format(
                    "BayesModel: genetic design {} rows {} != phenotype rows "
                    "{}",
                    genetic.type,
                    genetic.X.rows(),
                    num_individuals_));
        }
        if (std::ranges::contains(seen_modes, genetic.type))
        {
            throw GelexException(
                fmt::format(
                    "BayesModel: duplicate genetic mode {}", genetic.type));
        }
        seen_modes.push_back(genetic.type);
    }
}

auto BayesModel::visit(infra::FieldVisitor& visitor) const -> void
{
    auto state_scope = visitor.scope("state");
    {
        auto fixed_scope = visitor.scope("fixed");
        fixed_.visit(visitor);
    }
    for (auto [i, random] : std::views::enumerate(random_))
    {
        auto random_scope = visitor.scope(fmt::format("random_{}", i));
        random.visit(visitor);
    }
}

}  // namespace gelex
