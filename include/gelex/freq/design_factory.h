// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_DESIGN_FACTORY_H_
#define GELEX_FREQ_DESIGN_FACTORY_H_

#include <Eigen/Core>
#include <span>
#include <string>
#include <vector>

#include "gelex/data/dataframe/key_type.h"
#include "gelex/freq/design.h"

namespace gelex
{

template <KeyType Key>
class DataFrame;

template <KeyType Key>
class DataFrameIndex;

}  // namespace gelex

namespace gelex::freq
{

[[nodiscard]] auto make_random_designs(const DataFrame<std::string>& frame)
    -> std::vector<RandomDesign>;

[[nodiscard]] auto make_quantitative_random_design(
    const DataFrame<std::string>& frame,
    std::string name) -> RandomDesign;

[[nodiscard]] auto make_grm_designs(
    std::span<const std::string> prefixes,
    const DataFrameIndex<std::string>& index) -> std::vector<RandomDesign>;

[[nodiscard]] auto make_interaction_design(
    std::string name,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::MatrixXd>& rhs) -> RandomDesign;

}  // namespace gelex::freq

#endif  // GELEX_FREQ_DESIGN_FACTORY_H_
