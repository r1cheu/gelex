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

#include "gelex/algo/infer/posterior_calculator.h"

#include <Eigen/Core>

namespace gelex::posterior::detail
{

void compute_pve(
    PosteriorSummary& summary,
    const Eigen::Ref<const Eigen::VectorXd>& mean_coeffs,
    double phenotype_var)
{
    const Eigen::Index n_params = mean_coeffs.size();
    if (n_params == 0 || phenotype_var <= 0.0)
    {
        return;
    }

    summary.mean = mean_coeffs.array().square() / phenotype_var;
}

}  // namespace gelex::posterior::detail
