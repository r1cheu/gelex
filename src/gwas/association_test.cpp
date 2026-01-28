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

#include "gelex/gwas/association_test.h"
#include <mkl_vml_functions.h>

#include <Eigen/Dense>

namespace gelex::gwas
{

namespace
{

void chi2_test(
    Eigen::Ref<Eigen::VectorXd> wald_stats,
    Eigen::Ref<Eigen::VectorXd> p_values)
{
    int n = static_cast<int>(wald_stats.size());
    vdSqrt(n, wald_stats.data(), p_values.data());
    cblas_dscal(n, -1.0, p_values.data(), 1);
    vdCdfNorm(n, p_values.data(), p_values.data());
    cblas_dscal(n, 2.0, p_values.data(), 1);
}
}  // namespace

void wald_test(AssocInput& input, AssocOutput& output)
{
    output.zt_v_inv_r = (input.Z.transpose() * input.V_inv_y);
    input.W = input.V_inv * input.Z;
    output.zt_v_inv_z = (input.Z.transpose() * input.W).diagonal();

    output.beta = (output.zt_v_inv_r.array() / output.zt_v_inv_z.array());
    output.se = (1.0 / output.zt_v_inv_z.array()).sqrt();
    output.stats = (output.beta.array() / output.se.array()).square();
    chi2_test(output.stats, output.p_value);
}
}  // namespace gelex::gwas
