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

#ifndef GELEX_BAYES_GENETIC_DETAIL_APPLY_FITTED_UPDATE_H_
#define GELEX_BAYES_GENETIC_DETAIL_APPLY_FITTED_UPDATE_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <span>
#include <type_traits>
#include <variant>

#include "gelex/bayes/genotype/operations.h"
#include "gelex/bayes/genotype/projection.h"

namespace gelex::detail
{

template <std::size_t ExtraCount>
    requires(ExtraCount != std::dynamic_extent)
inline auto apply_fitted_update(
    const bayes::GeneticProjection& projection,
    Eigen::Index marker,
    const std::variant<
        std::monostate,
        bayes::AxpyTarget,
        std::array<bayes::AxpyTarget, 2>>& component_update,
    std::span<const bayes::AxpyTarget, ExtraCount> extra_targets) -> void
{
    std::array<bayes::AxpyTarget, ExtraCount + 2> targets{};
    std::size_t target_count = 0;
    for (const auto& update : extra_targets)
    {
        if (update.scale != 0.0)
        {
            targets[target_count++] = update;
        }
    }
    std::visit(
        [&](const auto& update)
        {
            using Update = std::remove_cvref_t<decltype(update)>;
            if constexpr (std::is_same_v<Update, bayes::AxpyTarget>)
                targets[target_count++] = update;
            else if constexpr (!std::is_same_v<Update, std::monostate>)
            {
                for (const auto& target : update)
                    targets[target_count++] = target;
            }
        },
        component_update);
    if (target_count != 0)
    {
        projection.axpy(
            marker,
            std::span<const bayes::AxpyTarget>{targets.data(), target_count});
    }
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_APPLY_FITTED_UPDATE_H_
