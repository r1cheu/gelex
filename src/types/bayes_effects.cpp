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

#include "../src/types/bayes_effects.h"

#include <utility>

#include <Eigen/Core>

namespace gelex
{
using Eigen::Ref;

using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

namespace bayes
{

FixedEffect::FixedEffect(
    std::optional<std::vector<std::string>> levels,
    MatrixXd&& design_matrix)
    : design_matrix(std::move(design_matrix)), levels(std::move(levels))
{
    cols_norm = this->design_matrix.colwise().squaredNorm();
}

FixedState::FixedState(const FixedEffect& effect)
    : coeffs(VectorXd::Zero(effect.design_matrix.cols())) {};

RandomEffect::RandomEffect(
    std::optional<std::vector<std::string>> levels,
    MatrixXd&& design_matrix)
    : design_matrix(std::move(design_matrix)), levels(std::move(levels))
{
    cols_norm = this->design_matrix.colwise().squaredNorm();
}

RandomState::RandomState(const RandomEffect& effect)
    : coeffs(VectorXd::Zero(effect.design_matrix.cols())),
      variance{effect.init_variance}
{
}

// BaseMarkerEffect constructors are inline in header

DominantState::DominantState(const DominantEffect& effect)
    : GeneticState{effect},
      ratios(VectorXd::Zero(get_cols(effect.design_matrix))),
      ratio_mean(effect.ratio_mean),
      ratio_variance(effect.ratio_variance)
{
    if (effect.init_pi)
    {
        tracker = VectorXi::Zero(get_cols(effect.design_matrix));
        pi
            = {effect.init_pi.value(),
               Eigen::VectorXi::Zero(effect.init_pi->size())};
    };
};

}  // namespace bayes
}  // namespace gelex
