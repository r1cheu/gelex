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

#include "gelex/freq/model.h"

#include <Eigen/Core>
#include <utility>
#include <vector>

#include "gelex/data/fixed_design.h"
#include "gelex/freq/design.h"
#include "gelex/infra/var.h"

namespace gelex
{

FreqModel::FreqModel(
    Eigen::VectorXd phenotype,
    FixedDesign fixed_design,
    std::vector<freq::RandomDesign> random)
    : phenotype_(std::move(phenotype)),
      phenotype_variance_(vecvar(phenotype_, VarNormType::Population)),
      fixed_(std::move(fixed_design)),
      random_(std::move(random))
{
    num_individuals_ = phenotype_.size();
}

FreqState::FreqState(const FreqModel& model)
    : Vp_(model.phenotype_variance()), fixed_(model.fixed())
{
    for (const auto& r : model.random())
    {
        random_.emplace_back(r);
    }
    init_variance_components(model);
}

auto FreqState::compute_variance_ratio() -> void
{
    double total_random_variance = 0.0;
    for (const auto& r : random_)
    {
        total_random_variance += r.variance;
    }

    Vp_ = total_random_variance + residual_.variance;

    if (Vp_ > 0.0)
    {
        for (auto& r : random_)
        {
            r.variance_ratio = r.variance / Vp_;
        }
    }
}

auto FreqState::init_variance_components(const FreqModel& model) -> void
{
    const double random_proportion = model.random().empty() ? 0.0 : 0.5;
    const double init_residual_variance = Vp_ * (1.0 - random_proportion);
    double init_random_variance = Vp_ * random_proportion;
    const auto num_random = static_cast<double>(model.random().size());

    init_random_variance
        = num_random < 1 ? 0.0 : init_random_variance / num_random;

    for (auto& r : random_)
    {
        r.variance = init_random_variance;
    }
    residual_.variance = init_residual_variance;
}

}  // namespace gelex
