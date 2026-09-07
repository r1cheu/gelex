// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0
#ifndef GELEX_DATA_RANK_INVERSE_NORM_TRANSFORM_H_
#define GELEX_DATA_RANK_INVERSE_NORM_TRANSFORM_H_

#include <Eigen/Dense>
#include <array>
#include <cstdint>
#include <string_view>
#include <utility>

namespace gelex
{
enum class RintType : std::uint8_t
{
    None,
    Direct,
    Indirect
};

inline constexpr std::array rint_type_names{
    std::pair{RintType::None, std::string_view{"none"}},
    std::pair{RintType::Direct, std::string_view{"dint"}},
    std::pair{RintType::Indirect, std::string_view{"iint"}},
};

auto direct_int(
    Eigen::Ref<Eigen::VectorXd> phenotype,
    double offset = 3.0 / 8.0) -> void;

auto indirect_int(
    Eigen::Ref<Eigen::VectorXd> phenotype,
    const Eigen::Ref<const Eigen::MatrixXd>& covariates,
    double offset = 3.0 / 8.0) -> void;

}  // namespace gelex

#endif  // GELEX_DATA_RANK_INVERSE_NORM_TRANSFORM_H_
