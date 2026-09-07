// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/data/fixed_design.h"

#include <Eigen/Core>
#include <Eigen/QR>
#include <cstddef>
#include <utility>

#include "gelex/data/dataframe/constants.h"
#include "gelex/exception.h"

namespace gelex
{

auto FixedDesign::make(Eigen::Index n_samples) -> FixedDesign
{
    return make(
        QuantitativeCovariate{.names = {}, .X = Eigen::MatrixXd(n_samples, 0)},
        DiscreteCovariate{.terms = {}, .X = Eigen::MatrixXd(n_samples, 0)});
}

auto FixedDesign::make(QuantitativeCovariate quantitative) -> FixedDesign
{
    const auto n_samples = quantitative.X.rows();
    return make(
        std::move(quantitative),
        DiscreteCovariate{.terms = {}, .X = Eigen::MatrixXd(n_samples, 0)});
}

auto FixedDesign::make(DiscreteCovariate discrete) -> FixedDesign
{
    const auto n_samples = discrete.X.rows();
    return make(
        QuantitativeCovariate{.names = {}, .X = Eigen::MatrixXd(n_samples, 0)},
        std::move(discrete));
}

auto FixedDesign::make(
    QuantitativeCovariate quantitative,
    DiscreteCovariate discrete) -> FixedDesign
{
    FixedDesign design;
    const auto n_samples = quantitative.X.rows();
    const auto quantitative_columns = quantitative.X.cols();
    const auto discrete_columns = discrete.X.cols();
    const auto total_columns
        = Eigen::Index{1} + quantitative_columns + discrete_columns;

    design.column_names_.reserve(static_cast<std::size_t>(total_columns));
    design.column_names_.emplace_back(intercept_name);
    design.column_names_.insert(
        design.column_names_.end(),
        quantitative.names.begin(),
        quantitative.names.end());
    for (const auto& term : discrete.terms)
    {
        for (const auto& level : term.levels)
        {
            if (level != term.reference_level)
            {
                design.column_names_.push_back(term.name + separator + level);
            }
        }
    }

    design.quantitative_names_ = std::move(quantitative.names);
    design.discrete_terms_ = std::move(discrete.terms);

    design.matrix_ = Eigen::MatrixXd::Zero(n_samples, total_columns);
    design.matrix_.col(0).setOnes();
    design.matrix_.middleCols(1, quantitative_columns) = quantitative.X;
    design.matrix_.middleCols(1 + quantitative_columns, discrete_columns)
        = discrete.X;
    design.xtx_diag_ = design.matrix_.colwise().squaredNorm();
    if (!design.matrix_.allFinite())
    {
        throw GelexException("FixedDesign: X must contain only finite values");
    }
    if (design.matrix_.fullPivHouseholderQr().rank() != total_columns)
    {
        throw GelexException("FixedDesign: X must have full column rank");
    }
    return design;
}

}  // namespace gelex
