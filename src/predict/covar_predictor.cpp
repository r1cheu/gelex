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

#include <Eigen/Core>
#include <cmath>
#include <format>

#include "gelex/exception.h"

namespace gelex
{

CovarPredictor::CovarPredictor(const detail::CovarEffects& effects)
    : effects_(&effects), data_(nullptr), n_samples_(0), n_cont_(0), n_cat_(0)
{
}

CovarPredictor::Result CovarPredictor::compute(const PredictData& data)
{
    data_ = &data;
    n_samples_ = data.genotype.rows();
    n_cont_ = data.qcovariate_names.size();
    n_cat_ = data.dcovariate_names.size();

    const size_t n_covariate = 1 + n_cont_ + n_cat_;

    result_.predictions = Eigen::MatrixXd::Zero(
        n_samples_, static_cast<Eigen::Index>(n_covariate));
    result_.names.clear();
    result_.names.reserve(n_covariate);

    validate_intercept();
    validate_qcovariates();
    validate_continuous_coefficients();
    validate_categorical_coefficients();

    compute_intercept();
    compute_continuous();
    compute_categorical();

    return result_;
}

void CovarPredictor::validate_intercept() const
{
    if (std::isnan(effects_->intercept))
    {
        throw InvalidInputException("Intercept coefficient is missing or NaN");
    }
}

void CovarPredictor::validate_qcovariates() const
{
    if (data_->qcovariates.cols() != static_cast<Eigen::Index>(n_cont_ + 1))
    {
        throw InvalidInputException(
            std::format(
                "qcovariates matrix has {} columns, expected {} ({} continuous "
                "+ intercept)",
                data_->qcovariates.cols(),
                n_cont_ + 1,
                n_cont_));
    }
}

void CovarPredictor::validate_continuous_coefficients() const
{
    for (size_t i = 0; i < n_cont_; ++i)
    {
        const std::string& var_name = data_->qcovariate_names[i];
        if (!effects_->continuous_coeffs.contains(var_name))
        {
            throw InvalidInputException(
                std::format(
                    "Missing coefficient for continuous variable '{}'",
                    var_name));
        }
    }
}

void CovarPredictor::validate_categorical_coefficients() const
{
    for (size_t i = 0; i < n_cat_; ++i)
    {
        const std::string& var_name = data_->dcovariate_names[i];
        auto cat_it = effects_->categorical_coeffs.find(var_name);
        if (cat_it == effects_->categorical_coeffs.end())
        {
            throw InvalidInputException(
                std::format(
                    "Missing coefficient for categorical variable '{}'",
                    var_name));
        }
    }
}

void CovarPredictor::compute_intercept()
{
    result_.predictions.col(0).setConstant(effects_->intercept);
    result_.names.emplace_back("Intercept");
}

void CovarPredictor::compute_continuous()
{
    for (size_t i = 0; i < n_cont_; ++i)
    {
        const auto col_i = static_cast<Eigen::Index>(1 + i);
        const std::string& var_name = data_->qcovariate_names[i];
        auto it = effects_->continuous_coeffs.find(var_name);
        result_.predictions.col(col_i)
            = data_->qcovariates.col(col_i) * it->second;
        result_.names.push_back(var_name);
    }
}

void CovarPredictor::compute_categorical()
{
    for (size_t i = 0; i < n_cat_; ++i)
    {
        const std::string& var_name = data_->dcovariate_names[i];
        auto cat_it = effects_->categorical_coeffs.find(var_name);

        const auto& level_coeffs = cat_it->second;
        const auto& levels_vec = data_->dcovariates.at(var_name);
        const auto cont_i = static_cast<Eigen::Index>(1 + n_cont_ + i);

        for (Eigen::Index j = 0; j < n_samples_; ++j)
        {
            const std::string& level = levels_vec[j];
            auto level_it = level_coeffs.find(level);
            if (level_it == level_coeffs.end())
            {
                throw InvalidInputException(
                    std::format(
                        "Missing coefficient for level '{}' of variable '{}'",
                        level,
                        var_name));
            }
            result_.predictions(j, cont_i) = level_it->second;
        }
        result_.names.push_back(var_name);
    }
}

}  // namespace gelex
