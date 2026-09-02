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
