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

#include "gelex/predict/compute.h"

#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include "gelex/io/predict_io.h"
#include "gelex/predict/types.h"

namespace gelex
{

auto compute_gebv(const GenotypeData& geno, const SnpEffects& effects)
    -> GEBVResult
{
    GEBVResult result;
    for (const auto& [mode, effect] : effects)
    {
        const auto geno_it = geno.find(mode);
        if (geno_it == geno.end())
        {
            continue;
        }
        Eigen::VectorXd component = geno_it->second * effect;
        if (result.total.size() == 0)
        {
            result.total = component;
        }
        else
        {
            result.total += component;
        }
        result.components.emplace(mode, std::move(component));
    }
    return result;
}

auto compute_covariate_effects(
    const Eigen::MatrixXd& covariates,
    const Coefficients& coefficients) -> CovariateResult
{
    Eigen::VectorXd total = covariates * coefficients.values;

    // Intercept is always the first entry; covariates start at index 1
    const auto n_samples = covariates.rows();
    const auto n_covars
        = static_cast<Eigen::Index>(coefficients.names.size()) - 1;

    std::vector<std::string> covar_names(
        coefficients.names.begin() + 1, coefficients.names.end());

    Eigen::MatrixXd per_covariate(n_samples, n_covars);
    for (Eigen::Index i = 0; i < n_covars; ++i)
    {
        per_covariate.col(i)
            = covariates.col(i + 1) * coefficients.values[i + 1];
    }

    return CovariateResult{
        .total = std::move(total),
        .per_covariate = std::move(per_covariate),
        .covar_names = std::move(covar_names)};
}

}  // namespace gelex
