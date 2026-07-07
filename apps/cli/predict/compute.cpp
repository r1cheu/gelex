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

#include "compute.h"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/data/dataframe/constants.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/encode.h"
#include "gelex/exception.h"

namespace cli
{

auto compute_gebv(
    const gelex::ModeMap<Eigen::MatrixXd>& geno,
    const gelex::ModeMap<Eigen::VectorXd>& effects)
    -> gelex::ModeMap<Eigen::VectorXd>
{
    gelex::ModeMap<Eigen::VectorXd> components;
    for (const auto& [mode, effect] : effects)
    {
        const auto geno_it = geno.find(mode);
        if (geno_it == geno.end())
        {
            continue;
        }
        components.emplace(mode, geno_it->second * effect);
    }
    return components;
}

auto build_covariate_design(
    std::span<const std::string> term_names,
    const std::optional<gelex::DataFrame<std::string>>& qcovar_df,
    const std::optional<gelex::DataFrame<std::string>>& dcovar_df,
    Eigen::Index n_samples) -> Eigen::MatrixXd
{
    // group dcovar terms by column name: "Sex\x1FM" -> col="Sex", level="M"
    std::map<std::string, std::vector<std::string>> dcovar_levels;
    for (const auto& name : term_names)
    {
        auto pos = name.find(gelex::SEPARATOR);
        if (pos != std::string::npos)
        {
            dcovar_levels[name.substr(0, pos)].push_back(name.substr(pos + 1));
        }
    }

    // encode dcovar columns and check levels against the fitted terms
    std::map<std::string, gelex::EncodedResult<>> encoded;
    for (const auto& [col_name, levels] : dcovar_levels)
    {
        if (!dcovar_df)
        {
            throw gelex::GelexException(
                fmt::format(
                    "discrete covariate file required for term '{}'",
                    col_name));
        }
        const auto& col = (*dcovar_df)[col_name];
        // TODO(rlchen): report level mismatch via observer
        [[maybe_unused]] auto mismatch = gelex::check_levels(col, levels);
        encoded[col_name] = gelex::encode(col, levels);
    }

    // build design matrix in term order
    const auto n_terms = static_cast<Eigen::Index>(term_names.size());
    Eigen::MatrixXd X(n_samples, n_terms);

    for (Eigen::Index i = 0; i < n_terms; ++i)
    {
        const auto& term = term_names[static_cast<std::size_t>(i)];

        if (term == "Intercept")
        {
            X.col(i).setOnes();
            continue;
        }

        auto sep_pos = term.find(gelex::SEPARATOR);
        if (sep_pos != std::string::npos)
        {
            auto col_name = term.substr(0, sep_pos);
            const auto& enc = encoded.at(col_name);
            auto it = std::ranges::find(enc.level_names, term);
            auto col_idx = std::distance(enc.level_names.begin(), it);
            X.col(i) = enc.data.col(col_idx);
        }
        else
        {
            if (!qcovar_df)
            {
                throw gelex::GelexException(
                    fmt::format(
                        "quantitative covariate file required for term '{}'",
                        term));
            }
            X.col(i) = (*qcovar_df)[term].to_map<double>();
        }
    }
    return X;
}

auto compute_covariate_effects(
    const Eigen::MatrixXd& covariates,
    std::span<const std::string> term_names,
    const Eigen::Ref<const Eigen::VectorXd>& coefficients) -> CovariateResult
{
    // Intercept is always the first entry; covariates start at index 1
    const auto n_samples = covariates.rows();
    const auto n_covars = static_cast<Eigen::Index>(term_names.size()) - 1;

    std::vector<std::string> covar_names(
        term_names.begin() + 1, term_names.end());

    Eigen::MatrixXd per_covariate(n_samples, n_covars);
    for (Eigen::Index i = 0; i < n_covars; ++i)
    {
        per_covariate.col(i) = covariates.col(i + 1) * coefficients[i + 1];
    }

    return CovariateResult{
        .per_covariate = std::move(per_covariate),
        .covar_names = std::move(covar_names)};
}

}  // namespace cli
