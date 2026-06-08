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

#include "gelex/data/covariates.h"

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/dataframe/encode.h"
#include "gelex/data/reader.h"

namespace gelex
{

auto make_quantitative_covariate(const dataframe::DataFrame<std::string>& frame)
    -> QuantitativeCovariate
{
    return QuantitativeCovariate{
        .names = frame.names(), .X = frame.to_mat<double>()};
}

auto make_discrete_covariate(const dataframe::DataFrame<std::string>& frame)
    -> DiscreteCovariate
{
    std::vector<std::string> names;
    std::vector<std::vector<std::string>> levels;
    std::vector<std::string> reference_levels;
    std::vector<dataframe::EncodedResult<>> encoded_results;

    for (std::size_t i = 0; i < frame.cols(); ++i)
    {
        const auto& col = frame.col(i);
        auto all_levels = dataframe::collect_levels(col);
        if (all_levels.size() < 2)
        {
            continue;
        }
        names.emplace_back(col.name());
        reference_levels.push_back(all_levels.front());
        encoded_results.push_back(
            dataframe::encode(
                col, std::span<const std::string>(all_levels).subspan(1)));
        levels.push_back(std::move(all_levels));
    }

    Eigen::Index total_cols = 0;
    for (const auto& r : encoded_results)
    {
        total_cols += r.data.cols();
    }

    Eigen::MatrixXd X(static_cast<Eigen::Index>(frame.rows()), total_cols);
    Eigen::Index col_offset = 0;
    for (const auto& r : encoded_results)
    {
        X.middleCols(col_offset, r.data.cols()) = r.data;
        col_offset += r.data.cols();
    }

    return DiscreteCovariate{
        .names = std::move(names),
        .levels = std::move(levels),
        .reference_levels = std::move(reference_levels),
        .X = std::move(X)};
}

auto make_random_designs(const dataframe::DataFrame<std::string>& frame)
    -> std::vector<freq::RandomDesign>
{
    std::vector<freq::RandomDesign> random_designs;
    random_designs.reserve(frame.cols());
    for (std::size_t i = 0; i < frame.cols(); ++i)
    {
        const auto& col = frame.col(i);
        auto result = dataframe::one_hot_encode(col);

        random_designs.emplace_back(
            result.name,
            std::move(result.level_names),
            std::nullopt,
            result.data * result.data.transpose());
    }
    return random_designs;
}

auto make_grm_designs(
    std::span<const std::string> prefixes,
    const dataframe::Index<std::string>& index)
    -> std::vector<freq::RandomDesign>
{
    std::vector<freq::RandomDesign> grm_designs;
    grm_designs.reserve(prefixes.size());
    for (const auto& prefix : prefixes)
    {
        auto term_name = prefix;
        auto K = read_grm(term_name, &index);
        grm_designs.emplace_back(
            term_name, std::nullopt, std::nullopt, std::move(K));
    }
    return grm_designs;
}
}  // namespace gelex
