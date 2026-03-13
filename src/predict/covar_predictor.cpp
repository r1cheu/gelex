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

#include "predict/covar_predictor.h"

#include <algorithm>
#include <format>

#include "gelex/exception.h"

namespace gelex
{

CovarPredictor::CovarPredictor(const CovariateParams& params) : params_(&params)
{
}

auto CovarPredictor::compute(const PredictData& data) -> Result
{
    const auto n_samples = data.genotype.rows();
    const auto n_terms = static_cast<Eigen::Index>(params_->term_names.size());

    Result result;
    result.predictions = Eigen::MatrixXd::Zero(n_samples, n_terms);
    result.names = params_->term_names;

    for (Eigen::Index i = 0; i < n_terms; ++i)
    {
        const auto& term = params_->term_names[static_cast<size_t>(i)];
        const double coeff = params_->coefficients(i);

        if (term == "Intercept")
        {
            result.predictions.col(i).setConstant(coeff);
            continue;
        }

        // Try continuous: find term in qcovariate_names
        auto cont_it = std::ranges::find(data.qcovariate_names, term);
        if (cont_it != data.qcovariate_names.end())
        {
            const auto col_idx = static_cast<Eigen::Index>(
                std::distance(data.qcovariate_names.begin(), cont_it));
            // qcovariates col 0 is intercept column, so offset by 1
            result.predictions.col(i)
                = data.qcovariates.col(col_idx + 1) * coeff;
            continue;
        }

        // Try categorical: parse "VarName_Level" pattern
        auto underscore_pos = term.find('_');
        if (underscore_pos != std::string::npos)
        {
            auto var_name = term.substr(0, underscore_pos);
            auto level = term.substr(underscore_pos + 1);

            auto dcovar_it = data.dcovariates.find(var_name);
            if (dcovar_it != data.dcovariates.end())
            {
                const auto& levels_vec = dcovar_it->second;
                for (Eigen::Index j = 0; j < n_samples; ++j)
                {
                    if (levels_vec[static_cast<size_t>(j)] == level)
                    {
                        result.predictions(j, i) = coeff;
                    }
                }
                continue;
            }
        }

        throw InvalidInputException(
            std::format("unknown term '{}' in covariate parameters", term));
    }

    return result;
}

}  // namespace gelex
