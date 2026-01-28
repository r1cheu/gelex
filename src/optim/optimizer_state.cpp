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

#include "gelex/optim/optimizer_state.h"

#include "gelex/model/freq/model.h"

namespace gelex
{

OptimizerState::OptimizerState(const FreqModel& model)
    : num_individuals_(model.num_individuals()),
      phenotype_variance_(model.phenotype_variance())
{
    v.resize(num_individuals_, num_individuals_);
    proj.resize(num_individuals_, num_individuals_);
    proj_y.resize(num_individuals_);
    tx_vinv_x.resize(model.fixed().X.cols(), model.fixed().X.cols());

    // preallocate for AI policy
    // n_comp = 1 (residual) + n_random + n_genetic
    auto n_comp = static_cast<Eigen::Index>(
        1 + model.random().size() + model.genetic().size());
    dvpy.resize(num_individuals_, n_comp);
    first_grad.resize(n_comp);
}

}  // namespace gelex
