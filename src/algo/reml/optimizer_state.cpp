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

#include "gelex/algo/reml/optimizer_state.h"
#include <Eigen/Core>

#include "gelex/freq/model.h"

namespace gelex::reml
{

OptimizerState::OptimizerState(const FreqModel& model)
    : num_individuals_(model.num_individuals()),
      phenotype_variance_(model.phenotype_variance())
{
    const auto n_fixed = model.fixed().X.cols();
    V.resize(num_individuals_, num_individuals_);
    Py.resize(num_individuals_);
    ViX.resize(num_individuals_, n_fixed);
    XtViX_inv.resize(n_fixed, n_fixed);

    // preallocate for AI policy
    // n_comp = 1 (residual) + n_random
    auto n_comp = static_cast<Eigen::Index>(1 + model.random().size());
    dvpy.resize(num_individuals_, n_comp);
    first_grad.resize(n_comp);
}

auto OptimizerState::trace_proj() const -> double
{
    return V.trace() - XtViX_inv.cwiseProduct(ViX.transpose() * ViX).sum();
}

auto OptimizerState::trace_proj_k(
    const Eigen::Ref<const Eigen::MatrixXd>& K) const -> double
{
    return V.cwiseProduct(K).sum()
           - XtViX_inv.cwiseProduct(ViX.transpose() * K * ViX).sum();
}

}  // namespace gelex::reml
