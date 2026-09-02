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

#ifndef GELEX_BAYES_GENETIC_COMPONENT_LAYOUT_H_
#define GELEX_BAYES_GENETIC_COMPONENT_LAYOUT_H_

#include <array>
#include <concepts>
#include <cstddef>

#include "gelex/types/genetic_mode.h"

namespace gelex
{

template <typename T>
concept ComponentLayout = requires(GeneticMode mode, std::size_t class_index) {
    { T::class_count } -> std::convertible_to<std::size_t>;
    { T::component_count } -> std::convertible_to<std::size_t>;
    { T::no_component } -> std::convertible_to<int>;
    { T::component_index(mode, class_index) } noexcept -> std::same_as<int>;
};

struct SingleComponentLayout
{
    static constexpr std::size_t class_count = 1;
    static constexpr std::size_t component_count = 1;
    static constexpr int no_component = -1;

    [[nodiscard]] static constexpr auto component_index(
        GeneticMode /*mode*/,
        std::size_t class_index) noexcept -> int
    {
        return class_index == 0 ? 0 : no_component;
    }
};

template <std::size_t ClassCount>
    requires(ClassCount > 1)
struct ZeroInflatedComponentLayout
{
    static constexpr std::size_t class_count = ClassCount;
    static constexpr std::size_t component_count = ClassCount - 1;
    static constexpr int no_component = -1;

    [[nodiscard]] static constexpr auto component_index(
        GeneticMode /*mode*/,
        std::size_t class_index) noexcept -> int
    {
        return class_index > 0 && class_index < class_count
                   ? static_cast<int>(class_index - 1)
                   : no_component;
    }
};

struct JointZeroInflatedComponentLayout
{
    static constexpr std::size_t class_count = 4;
    static constexpr std::size_t component_count = 2;
    static constexpr int no_component = -1;

    [[nodiscard]] static constexpr auto component_index(
        GeneticMode mode,
        std::size_t class_index) noexcept -> int
    {
        if (class_index >= class_count)
        {
            return no_component;
        }
        switch (mode)
        {
            case GeneticMode::A:
                return additive_component_indices.at(class_index);
            case GeneticMode::D:
                return dominance_component_indices.at(class_index);
        }
        return no_component;
    }

   private:
    static constexpr std::array<int, class_count> additive_component_indices{
        no_component,
        0,
        no_component,
        1};
    static constexpr std::array<int, class_count> dominance_component_indices{
        no_component,
        no_component,
        0,
        1};
};

static_assert(ComponentLayout<SingleComponentLayout>);
static_assert(ComponentLayout<ZeroInflatedComponentLayout<2>>);
static_assert(ComponentLayout<JointZeroInflatedComponentLayout>);

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_COMPONENT_LAYOUT_H_
