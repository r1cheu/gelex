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

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/state.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class DenseResidualAdjustmentGuard
{
   public:
    DenseResidualAdjustmentGuard(
        Eigen::Ref<const Eigen::VectorXd> column,
        double& coeff,
        bayes::ResidualState& residual)
        : column_(std::move(column)),
          coeff_(coeff),
          residual_(residual),
          old_value_(coeff)
    {
    }

    DenseResidualAdjustmentGuard(const DenseResidualAdjustmentGuard&) = delete;
    DenseResidualAdjustmentGuard(DenseResidualAdjustmentGuard&&) = delete;
    auto operator=(const DenseResidualAdjustmentGuard&)
        -> DenseResidualAdjustmentGuard& = delete;
    auto operator=(DenseResidualAdjustmentGuard&&)
        -> DenseResidualAdjustmentGuard& = delete;

    ~DenseResidualAdjustmentGuard()
    {
        const double diff = old_value_ - coeff_;
        if (diff != 0.0)
        {
            residual_.y_adj.array() += diff * column_.array();
        }
    }

   private:
    Eigen::Ref<const Eigen::VectorXd> column_;
    double& coeff_;
    bayes::ResidualState& residual_;
    double old_value_{0.0};
};

class GeneticResidualAdjustmentGuard
{
   public:
    GeneticResidualAdjustmentGuard(
        const bayes::GeneticDesign& design,
        GeneticMode mode,
        Eigen::Index marker,
        bayes::GeneticState& state,
        bayes::ResidualState& residual)
        : design_(design),
          mode_(mode),
          marker_(marker),
          state_(state),
          residual_(residual),
          old_value_(state.coeffs(marker))
    {
    }

    GeneticResidualAdjustmentGuard(const GeneticResidualAdjustmentGuard&)
        = delete;
    GeneticResidualAdjustmentGuard(GeneticResidualAdjustmentGuard&&) = delete;
    auto operator=(const GeneticResidualAdjustmentGuard&)
        -> GeneticResidualAdjustmentGuard& = delete;
    auto operator=(GeneticResidualAdjustmentGuard&&)
        -> GeneticResidualAdjustmentGuard& = delete;

    ~GeneticResidualAdjustmentGuard()
    {
        const double diff = old_value_ - state_.coeffs(marker_);
        if (diff != 0.0)
        {
            design_.axpy(mode_, marker_, diff, residual_.y_adj);
        }
    }

   private:
    const bayes::GeneticDesign& design_;
    GeneticMode mode_;
    Eigen::Index marker_;
    bayes::GeneticState& state_;
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
        const bayes::GeneticDesign& design,
        GeneticMode mode,
        Eigen::Index marker,
        bayes::GeneticState& state,
        bayes::SingleScaledMixtureGaussianState& prior_state)
        : design_(design),
          mode_(mode),
          marker_(marker),
          state_(state),
          prior_state_(prior_state),
          old_value_(state.coeffs(marker)),
          old_class_(prior_state.assignment()(marker))
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
        const double coeff = state_.coeffs(marker_);
        auto& component = prior_state_.component();
        const int old_class = old_class_;
        const int new_class = prior_state_.assignment()(marker_);
        if (old_class == new_class)
        {
            if (old_class > 0 && old_value_ != coeff)
            {
                design_.axpy(
                    mode_,
                    marker_,
                    coeff - old_value_,
                    component.gebv[old_class - 1]);
            }
            return;
        }
        if (old_class > 0)
        {
            design_.axpy(
                mode_, marker_, -old_value_, component.gebv[old_class - 1]);
        }
        if (new_class > 0)
        {
            design_.axpy(mode_, marker_, coeff, component.gebv[new_class - 1]);
        }
    }

   private:
    const bayes::GeneticDesign& design_;
    GeneticMode mode_;
    Eigen::Index marker_;
    bayes::GeneticState& state_;
    bayes::SingleScaledMixtureGaussianState& prior_state_;
    double old_value_{0.0};
    int old_class_{0};
};

class JointGeneticAdjustmentGuard
{
   public:
    JointGeneticAdjustmentGuard(
        const bayes::GeneticDesign& design,
        Eigen::Index marker,
        bayes::JointGeneticBlockState& block,
        bayes::ResidualState& residual)
        : design_(design),
          marker_(marker),
          block_(block),
          residual_(residual),
          old_additive_(block.state(GeneticMode::A).coeffs(marker)),
          old_dominance_(block.state(GeneticMode::D).coeffs(marker))
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
        auto& additive = block_.state(GeneticMode::A);
        auto& dominance = block_.state(GeneticMode::D);
        const double additive_diff = old_additive_ - additive.coeffs(marker_);
        const double dominance_diff
            = old_dominance_ - dominance.coeffs(marker_);
        if (additive_diff == 0.0 && dominance_diff == 0.0)
        {
            return;
        }
        if (additive_diff != 0.0)
        {
            design_.axpy(
                GeneticMode::A, marker_, additive_diff, residual_.y_adj);
            design_.axpy(GeneticMode::A, marker_, -additive_diff, additive.u);
        }
        if (dominance_diff != 0.0)
        {
            design_.axpy(
                GeneticMode::D, marker_, dominance_diff, residual_.y_adj);
            design_.axpy(GeneticMode::D, marker_, -dominance_diff, dominance.u);
        }
    }

   private:
    const bayes::GeneticDesign& design_;
    Eigen::Index marker_;
    bayes::JointGeneticBlockState& block_;
    bayes::ResidualState& residual_;
    double old_additive_{0.0};
    double old_dominance_{0.0};
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex

#endif  // GELEX_ALGO_MCMC_INVARIANT_H_
