// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/data/covariates.h"

#include <Eigen/Core>
#include <cstddef>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/encode.h"

namespace gelex
{

auto make_quantitative_covariate(const DataFrame<std::string>& frame)
    -> QuantitativeCovariate
{
    return QuantitativeCovariate{
        .names = frame.names(), .X = frame.to_mat<double>()};
}

auto make_discrete_covariate(const DataFrame<std::string>& frame)
    -> DiscreteCovariate
{
    std::vector<DiscreteCovariateTerm> terms;
    std::vector<EncodedResult<>> encoded_results;

    for (std::size_t i = 0; i < frame.cols(); ++i)
    {
        const auto& col = frame.col(i);
        auto all_levels = collect_levels(col);
        if (all_levels.size() < 2)
        {
            continue;
        }
        auto reference_level = all_levels.front();
        encoded_results.push_back(
            encode(col, std::span<const std::string>(all_levels).subspan(1)));
        terms.push_back(
            DiscreteCovariateTerm{
                .name = std::string(col.name()),
                .levels = std::move(all_levels),
                .reference_level = std::move(reference_level)});
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

    return DiscreteCovariate{.terms = std::move(terms), .X = std::move(X)};
}
}  // namespace gelex
