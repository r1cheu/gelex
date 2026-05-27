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

#ifndef GELEX_ALGO_INFER_MCMC_INVARIANT_H_
#define GELEX_ALGO_INFER_MCMC_INVARIANT_H_

#include <Eigen/Core>

#include "gelex/model/bayes/state.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class ResidualAdjustmentGuard
{
   public:
    ResidualAdjustmentGuard(
        Eigen::Ref<const Eigen::VectorXd> column,
        double& coeff,
        bayes::ResidualState& residual)
        : column_(column), coeff_(coeff), residual_(residual), old_value_(coeff)
    {
    }

    ResidualAdjustmentGuard(const ResidualAdjustmentGuard&) = delete;
    ResidualAdjustmentGuard(ResidualAdjustmentGuard&&) = delete;
    auto operator=(const ResidualAdjustmentGuard&)
        -> ResidualAdjustmentGuard& = delete;
    auto operator=(ResidualAdjustmentGuard&&)
        -> ResidualAdjustmentGuard& = delete;

    ~ResidualAdjustmentGuard()
    {
        const double diff = old_value_ - coeff_;
        if (diff == 0.0)
        {
            return;
        }
        residual_.y_adj.array() += diff * column_.array();
    }

   private:
    Eigen::Ref<const Eigen::VectorXd> column_;
    double& coeff_;
    bayes::ResidualState& residual_;
    double old_value_{0.0};
};

class GeneticAdjustmentGuard
{
   public:
    GeneticAdjustmentGuard(
        Eigen::Ref<const Eigen::VectorXd> column,
        double& coeff,
        bayes::ResidualState& residual,
        bayes::GeneticState& state)
        : column_(column),
          coeff_(coeff),
          residual_(residual),
          state_(state),
          old_value_(coeff)
    {
    }

    GeneticAdjustmentGuard(const GeneticAdjustmentGuard&) = delete;
    GeneticAdjustmentGuard(GeneticAdjustmentGuard&&) = delete;
    auto operator=(const GeneticAdjustmentGuard&)
        -> GeneticAdjustmentGuard& = delete;
    auto operator=(GeneticAdjustmentGuard&&)
        -> GeneticAdjustmentGuard& = delete;

    ~GeneticAdjustmentGuard()
    {
        const double diff = old_value_ - coeff_;
        if (diff == 0.0)
        {
            return;
        }
        for (Eigen::Index i = 0; i < column_.size(); ++i)
        {
            const double delta = diff * column_(i);
            residual_.y_adj(i) += delta;
            state_.u(i) -= delta;
        }
    }

   private:
    Eigen::Ref<const Eigen::VectorXd> column_;
    double& coeff_;
    bayes::ResidualState& residual_;
    bayes::GeneticState& state_;
    double old_value_{0.0};
};

class GeneticMixtureAdjustmentGuard
{
   public:
    GeneticMixtureAdjustmentGuard(
        Eigen::Ref<const Eigen::VectorXd> column,
        double& coeff,
        bayes::ResidualState& residual,
        bayes::GeneticState& state,
        bayes::ComponentState& component,
        Eigen::VectorXi& assignment,
        Eigen::Index marker_index)
        : column_(column),
          coeff_(coeff),
          residual_(residual),
          state_(state),
          component_(component),
          assignment_(assignment),
          marker_index_(marker_index),
          old_value_(coeff),
          old_class_(assignment(marker_index))
    {
    }

    GeneticMixtureAdjustmentGuard(const GeneticMixtureAdjustmentGuard&)
        = delete;
    GeneticMixtureAdjustmentGuard(GeneticMixtureAdjustmentGuard&&) = delete;
    auto operator=(const GeneticMixtureAdjustmentGuard&)
        -> GeneticMixtureAdjustmentGuard& = delete;
    auto operator=(GeneticMixtureAdjustmentGuard&&)
        -> GeneticMixtureAdjustmentGuard& = delete;

    ~GeneticMixtureAdjustmentGuard()
    {
        const double diff = old_value_ - coeff_;
        const int old_class = old_class_;
        const int new_class = assignment_(marker_index_);
        if (diff == 0.0 && old_class == new_class)
        {
            return;
        }
        for (Eigen::Index i = 0; i < column_.size(); ++i)
        {
            const double xi = column_(i);
            const double delta = diff * xi;
            residual_.y_adj(i) += delta;
            state_.u(i) -= delta;
            if (old_class > 0)
            {
                component_.gebv[old_class - 1](i) -= old_value_ * xi;
            }
            if (new_class > 0)
            {
                component_.gebv[new_class - 1](i) += coeff_ * xi;
            }
        }
    }

   private:
    Eigen::Ref<const Eigen::VectorXd> column_;
    double& coeff_;
    bayes::ResidualState& residual_;
    bayes::GeneticState& state_;
    bayes::ComponentState& component_;
    Eigen::VectorXi& assignment_;
    Eigen::Index marker_index_{};
    double old_value_{0.0};
    int old_class_{0};
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_INVARIANT_H_
