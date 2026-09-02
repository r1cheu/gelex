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

#include "gelex/freq/design_factory.h"

#include <Eigen/Core>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/encode.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/grm/io.h"
#include "gelex/freq/design.h"

namespace gelex::freq
{

auto make_random_designs(const DataFrame<std::string>& frame)
    -> std::vector<RandomDesign>
{
    std::vector<RandomDesign> random_designs;
    random_designs.reserve(frame.cols());
    for (std::size_t i = 0; i < frame.cols(); ++i)
    {
        const auto& col = frame.col(i);
        auto result = one_hot_encode(col);

        random_designs.emplace_back(
            result.name,
            std::move(result.level_names),
            std::nullopt,
            result.data * result.data.transpose(),
            RandomKind::Discrete);
    }
    return random_designs;
}

auto make_quantitative_random_design(
    const DataFrame<std::string>& frame,
    std::string name) -> RandomDesign
{
    Eigen::MatrixXd Z = frame.to_mat<double>();
    return RandomDesign{
        .name = std::move(name),
        .levels = std::nullopt,
        .Z = std::nullopt,
        .K = Z * Z.transpose(),
        .kind = RandomKind::Quantitative};
}

auto make_grm_designs(
    std::span<const std::string> prefixes,
    const DataFrameIndex<std::string>& index) -> std::vector<RandomDesign>
{
    std::vector<RandomDesign> grm_designs;
    grm_designs.reserve(prefixes.size());
    for (const auto& prefix : prefixes)
    {
        auto name = std::filesystem::path(prefix).filename().string();
        auto K = read_grm(prefix, &index);
        grm_designs.emplace_back(
            name, std::nullopt, std::nullopt, std::move(K), RandomKind::Grm);
    }
    return grm_designs;
}

auto make_interaction_design(
    std::string name,
    const Eigen::Ref<const Eigen::MatrixXd>& lhs,
    const Eigen::Ref<const Eigen::MatrixXd>& rhs) -> RandomDesign
{
    Eigen::MatrixXd K = lhs.cwiseProduct(rhs);
    K /= K.diagonal().mean();
    return RandomDesign{
        .name = std::move(name),
        .levels = std::nullopt,
        .Z = std::nullopt,
        .K = std::move(K),
        .kind = RandomKind::Interaction};
}

}  // namespace gelex::freq
