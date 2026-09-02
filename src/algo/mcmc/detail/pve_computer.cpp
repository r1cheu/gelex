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

#include "algo/mcmc/detail/pve_computer.h"

#include <Eigen/Core>

#include "gelex/bayes/design.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

PveComputer::PveComputer(
    const bayes::GeneticDesign& design,
    double phenotype_var)
    : design_(design), phenotype_var_(phenotype_var)
{
    if (phenotype_var_ <= 0.0)
    {
        throw GelexException(
            "PveComputer: phenotype variance must be positive");
    }

    if (!design_.modes().contains(GeneticMode::A)
        || !design_.modes().contains(GeneticMode::D))
    {
        return;
    }

    cov_ad_ = design_.col_covariance(GeneticMode::A, GeneticMode::D);
}

auto PveComputer::single(
    GeneticMode mode,
    const Eigen::Ref<const Eigen::VectorXd>& beta) const -> Eigen::VectorXd
{
    const auto& col_var = design_.col_var(mode);
    if (col_var.size() != beta.size())
    {
        throw GelexException("PveComputer: beta and col_var size mismatch");
    }

    return (col_var.transpose().array() * beta.array().square()
            / phenotype_var_)
        .matrix()
        .eval();
}

auto PveComputer::total(
    const Eigen::Ref<const Eigen::VectorXd>& beta_A,
    const Eigen::Ref<const Eigen::VectorXd>& beta_D) const -> Eigen::VectorXd
{
    const auto& additive_col_var = design_.col_var(GeneticMode::A);
    const auto& dominance_col_var = design_.col_var(GeneticMode::D);
    if (beta_A.size() != beta_D.size()
        || beta_A.size() != additive_col_var.size()
        || beta_A.size() != dominance_col_var.size()
        || beta_A.size() != cov_ad_.size())
    {
        throw GelexException("PveComputer: total A/D size mismatch");
    }

    return (beta_A.array().square() * additive_col_var.transpose().array()
            + beta_D.array().square() * dominance_col_var.transpose().array()
            + 2.0 * beta_A.array() * beta_D.array()
                  * cov_ad_.transpose().array())
               .matrix()
           / phenotype_var_;
}

}  // namespace gelex::detail
