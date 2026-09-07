// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_COVARIATES_H_
#define GELEX_DATA_COVARIATES_H_

#include <Eigen/Core>
#include <string>
#include <vector>

#include "gelex/data/dataframe/key_type.h"

namespace gelex
{

template <KeyType Key>
class DataFrame;

struct QuantitativeCovariate
{
    std::vector<std::string> names;
    Eigen::MatrixXd X;
};

struct DiscreteCovariateTerm
{
    std::string name;
    std::vector<std::string> levels;
    std::string reference_level;
};

struct DiscreteCovariate
{
    std::vector<DiscreteCovariateTerm> terms;
    Eigen::MatrixXd X;
};

[[nodiscard]] auto make_quantitative_covariate(
    const DataFrame<std::string>& frame) -> QuantitativeCovariate;

[[nodiscard]] auto make_discrete_covariate(const DataFrame<std::string>& frame)
    -> DiscreteCovariate;

}  // namespace gelex

#endif  // GELEX_DATA_COVARIATES_H_
