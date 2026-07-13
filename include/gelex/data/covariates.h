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
#include <span>
#include <string>
#include <vector>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/freq/design.h"
#include "gelex/types/fixed_designs.h"

namespace gelex
{

[[nodiscard]] auto make_quantitative_covariate(
    const DataFrame<std::string>& frame) -> QuantitativeCovariate;

[[nodiscard]] auto make_discrete_covariate(const DataFrame<std::string>& frame)
    -> DiscreteCovariate;

[[nodiscard]] auto make_random_designs(const DataFrame<std::string>& frame)
    -> std::vector<freq::RandomDesign>;

[[nodiscard]] auto make_quantitative_random_design(
    const DataFrame<std::string>& frame,
    std::string name) -> freq::RandomDesign;

[[nodiscard]] auto make_grm_designs(
    std::span<const std::string> prefixes,
    const DataFrameIndex<std::string>&) -> std::vector<freq::RandomDesign>;

// Hadamard product of two aligned kernels, rescaled to unit mean diagonal so
// the interaction variance component shares the scale of its base components.
[[nodiscard]] auto make_interaction_design(
    std::string name,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::MatrixXd>& rhs) -> freq::RandomDesign;

}  // namespace gelex

#endif  // GELEX_DATA_COVARIATES_H_
