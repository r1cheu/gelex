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
#ifndef GELEX_INFRA_STATS_RANK_INVERSE_NORM_TRANSFORM_H_
#define GELEX_INFRA_STATS_RANK_INVERSE_NORM_TRANSFORM_H_

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

inline constexpr std::array RINT_TYPE_NAMES{
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

#endif  // GELEX_INFRA_STATS_RANK_INVERSE_NORM_TRANSFORM_H_
