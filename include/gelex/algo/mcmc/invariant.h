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

#ifndef GELEX_ALGO_MCMC_INVARIANT_H_
#define GELEX_ALGO_MCMC_INVARIANT_H_

#include <Eigen/Core>
#include <utility>

#include "gelex/bayes/state.h"

namespace gelex
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class ResidualAdjustmentGuard
{
   public:
    ResidualAdjustmentGuard(
        Eigen::Ref<const Eigen::VectorXd> column,
        double& coeff,
        bayes::ResidualState& residual)
        : column_(std::move(column)),
          coeff_(coeff),
          residual_(residual),
          old_value_(coeff)
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

class GeneticSweepAdjustmentGuard
{
   public:
    GeneticSweepAdjustmentGuard(
        bayes::ResidualState& residual,
        bayes::GeneticState& state)
        : residual_(residual), state_(state)
    {
        residual_.old_y_adj = residual_.y_adj;
    }

    GeneticSweepAdjustmentGuard(const GeneticSweepAdjustmentGuard&) = delete;
    GeneticSweepAdjustmentGuard(GeneticSweepAdjustmentGuard&&) = delete;
    auto operator=(const GeneticSweepAdjustmentGuard&)
        -> GeneticSweepAdjustmentGuard& = delete;
    auto operator=(GeneticSweepAdjustmentGuard&&)
        -> GeneticSweepAdjustmentGuard& = delete;

    ~GeneticSweepAdjustmentGuard()
    {
        state_.u.noalias() += residual_.old_y_adj - residual_.y_adj;
    }

   private:
    bayes::ResidualState& residual_;
    bayes::GeneticState& state_;
};

class ComponentGebvAdjustmentGuard
{
   public:
    ComponentGebvAdjustmentGuard(
        Eigen::Ref<const Eigen::VectorXd> column,
        double& coeff,
        bayes::ComponentState& component,
        int& assignment)
        : column_(std::move(column)),
          coeff_(coeff),
          component_(component),
          assignment_(assignment),
          old_value_(coeff),
          old_class_(assignment)
    {
    }

    ComponentGebvAdjustmentGuard(const ComponentGebvAdjustmentGuard&) = delete;
    ComponentGebvAdjustmentGuard(ComponentGebvAdjustmentGuard&&) = delete;
    auto operator=(const ComponentGebvAdjustmentGuard&)
        -> ComponentGebvAdjustmentGuard& = delete;
    auto operator=(ComponentGebvAdjustmentGuard&&)
        -> ComponentGebvAdjustmentGuard& = delete;

    ~ComponentGebvAdjustmentGuard()
    {
        const int old_class = old_class_;
        const int new_class = assignment_;
        if (old_class == new_class)
        {
            if (old_class > 0 && old_value_ != coeff_)
            {
                component_.gebv[old_class - 1].noalias()
                    += (coeff_ - old_value_) * column_;
            }
            return;
        }
        if (old_class > 0)
        {
            component_.gebv[old_class - 1].noalias() += (-old_value_) * column_;
        }
        if (new_class > 0)
        {
            component_.gebv[new_class - 1].noalias() += coeff_ * column_;
        }
    }

   private:
    Eigen::Ref<const Eigen::VectorXd> column_;
    double& coeff_;
    bayes::ComponentState& component_;
    int& assignment_;
    double old_value_{0.0};
    int old_class_{0};
};

class JointGeneticAdjustmentGuard
{
   public:
    JointGeneticAdjustmentGuard(
        Eigen::Ref<const Eigen::VectorXd> first_column,
        Eigen::Ref<const Eigen::VectorXd> second_column,
        double& first_coeff,
        double& second_coeff,
        bayes::ResidualState& residual,
        bayes::GeneticState& first_state,
        bayes::GeneticState& second_state)
        : first_column_(std::move(first_column)),
          second_column_(std::move(second_column)),
          first_coeff_(first_coeff),
          second_coeff_(second_coeff),
          residual_(residual),
          first_state_(first_state),
          second_state_(second_state),
          old_first_(first_coeff),
          old_second_(second_coeff)
    {
    }

    JointGeneticAdjustmentGuard(const JointGeneticAdjustmentGuard&) = delete;
    JointGeneticAdjustmentGuard(JointGeneticAdjustmentGuard&&) = delete;
    auto operator=(const JointGeneticAdjustmentGuard&)
        -> JointGeneticAdjustmentGuard& = delete;
    auto operator=(JointGeneticAdjustmentGuard&&)
        -> JointGeneticAdjustmentGuard& = delete;

    ~JointGeneticAdjustmentGuard()
    {
        const double first_diff = old_first_ - first_coeff_;
        const double second_diff = old_second_ - second_coeff_;
        if (first_diff == 0.0 && second_diff == 0.0)
        {
            return;
        }
        if (first_diff != 0.0)
        {
            residual_.y_adj.noalias() += first_diff * first_column_;
            first_state_.u.noalias() += (-first_diff) * first_column_;
        }
        if (second_diff != 0.0)
        {
            residual_.y_adj.noalias() += second_diff * second_column_;
            second_state_.u.noalias() += (-second_diff) * second_column_;
        }
    }

   private:
    Eigen::Ref<const Eigen::VectorXd> first_column_;
    Eigen::Ref<const Eigen::VectorXd> second_column_;
    double& first_coeff_;
    double& second_coeff_;
    bayes::ResidualState& residual_;
    bayes::GeneticState& first_state_;
    bayes::GeneticState& second_state_;
    double old_first_{0.0};
    double old_second_{0.0};
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex

#endif  // GELEX_ALGO_MCMC_INVARIANT_H_
