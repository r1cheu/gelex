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

#include "gelex/bayes/model.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::mcmc::detail
{

PveComputer::PveComputer(const BayesModel& model, double phenotype_var)
    : model_(model), phenotype_var_(phenotype_var)
{
    if (phenotype_var_ <= 0.0)
    {
        throw GelexException(
            "PveComputer: phenotype variance must be positive");
    }

    const auto* additive = model_.genetic(GeneticMode::A);
    const auto* dominance = model_.genetic(GeneticMode::D);
    if (additive == nullptr || dominance == nullptr)
    {
        return;
    }

    if (additive->X.cols() != dominance->X.cols())
    {
        throw GelexException("PveComputer: joint A/D design size mismatch");
    }

    cov_ad_ = Eigen::RowVectorXd::Zero(additive->X.cols());
    for (Eigen::Index i = 0; i < additive->X.cols(); ++i)
    {
        cov_ad_(i)
            = additive->X.matrix().col(i).dot(dominance->X.matrix().col(i))
              / static_cast<double>(additive->X.rows());
    }
}

auto PveComputer::single(
    GeneticMode mode,
    const Eigen::Ref<const Eigen::VectorXd>& beta) const -> Eigen::VectorXd
{
    const auto* design = model_.genetic(mode);
    if (design == nullptr)
    {
        throw GelexException("PveComputer: missing genetic design");
    }

    if (design->col_var.size() != beta.size())
    {
        throw GelexException("PveComputer: beta and col_var size mismatch");
    }

    return (design->col_var.transpose().array() * beta.array().square()
            / phenotype_var_)
        .matrix()
        .eval();
}

auto PveComputer::total(
    const Eigen::Ref<const Eigen::VectorXd>& beta_A,
    const Eigen::Ref<const Eigen::VectorXd>& beta_D) const -> Eigen::VectorXd
{
    const auto* additive = model_.genetic(GeneticMode::A);
    const auto* dominance = model_.genetic(GeneticMode::D);
    if (additive == nullptr || dominance == nullptr)
    {
        throw GelexException("PveComputer: missing A/D genetic design");
    }

    if (beta_A.size() != beta_D.size()
        || beta_A.size() != additive->col_var.size()
        || beta_A.size() != dominance->col_var.size()
        || beta_A.size() != cov_ad_.size())
    {
        throw GelexException("PveComputer: total A/D size mismatch");
    }

    return (beta_A.array().square() * additive->col_var.transpose().array()
            + beta_D.array().square() * dominance->col_var.transpose().array()
            + 2.0 * beta_A.array() * beta_D.array()
                  * cov_ad_.transpose().array())
               .matrix()
           / phenotype_var_;
}

}  // namespace gelex::mcmc::detail
