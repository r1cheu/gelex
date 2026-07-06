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

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include "gelex/io/predict_reader.h"
#include "gelex/predict/types.h"

namespace gelex
{

auto compute_gebv(const GenotypeData& geno, const SnpEffects& effects)
    -> GEBVResult
{
    std::optional<Eigen::VectorXd> add_predictions;
    if (geno.add && effects.add)
    {
        add_predictions = (*geno.add) * (*effects.add);
    }

    std::optional<Eigen::VectorXd> dom_predictions;
    if (geno.dom && effects.dom)
    {
        dom_predictions = (*geno.dom) * (*effects.dom);
    }

    Eigen::VectorXd total;
    if (add_predictions && dom_predictions)
    {
        total = *add_predictions + *dom_predictions;
    }
    else if (add_predictions)
    {
        total = *add_predictions;
    }
    else if (dom_predictions)
    {
        total = *dom_predictions;
    }

    return GEBVResult{
        .total = std::move(total),
        .add_predictions = std::move(add_predictions),
        .dom_predictions = std::move(dom_predictions)};
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
