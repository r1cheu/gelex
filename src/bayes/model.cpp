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

#include <Eigen/Core>
#include <cmath>
#include <fmt/format.h>
#include <utility>
#include <vector>

#include "gelex/bayes/design_data.h"
#include "gelex/exception.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/types/fixed_designs.h"

namespace gelex
{
BayesModel::BayesModel(
    Eigen::VectorXd phenotype,
    FixedDesign fixed_design,
    std::vector<bayes::RandomDesign> random,
    bayes::GeneticDesign genetic)
    : phenotype_(std::move(phenotype)),
      fixed_(std::move(fixed_design)),
      random_(std::move(random)),
      genetic_(std::move(genetic))
{
    num_individuals_ = phenotype_.rows();
    if (num_individuals_ == 0)
    {
        throw GelexException("BayesModel: phenotype must not be empty");
    }

    phenotype_var_
        = detail::vecvar(phenotype_, detail::VarNormType::Population);
    if (!std::isfinite(phenotype_var_) || phenotype_var_ <= 0.0)
    {
        throw GelexException(
            fmt::format(
                "BayesModel: phenotype variance must be finite and positive, "
                "got {}",
                phenotype_var_));
    }

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

    if (genetic_.rows() != num_individuals_)
    {
        throw GelexException(
            fmt::format(
                "BayesModel: genetic design rows {} != phenotype rows {}",
                genetic_.rows(),
                num_individuals_));
    }
}

}  // namespace gelex
