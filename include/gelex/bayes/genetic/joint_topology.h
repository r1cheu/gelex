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

#ifndef GELEX_BAYES_GENETIC_JOINT_TOPOLOGY_H_
#define GELEX_BAYES_GENETIC_JOINT_TOPOLOGY_H_

#include <type_traits>
#include <utility>

#include "gelex/bayes/genetic/independent_topology.h"

namespace gelex
{

namespace detail
{

template <typename T>
inline constexpr bool is_independent_topology = false;

template <GeneticModeSet Modes, typename... Ts>
inline constexpr bool is_independent_topology<IndependentTopology<Modes, Ts...>>
    = true;

}  // namespace detail

template <typename ModeTopology, typename JointT>
    requires(
        detail::is_independent_topology<ModeTopology>
        && ModeTopology::modes == (GeneticMode::A | GeneticMode::D)
        && std::is_object_v<JointT>)
class JointTopology
{
   public:
    using mode_topology_type = ModeTopology;
    using joint_value_type = JointT;

    template <GeneticMode Mode>
    using mode_value_type =
        typename mode_topology_type::template mode_value_type<Mode>;

    static constexpr GeneticModeSet modes = mode_topology_type::modes;

    constexpr JointTopology() = default;

    constexpr JointTopology(mode_topology_type mode_values, JointT joint)
        : mode_values_{std::move(mode_values)}, joint_{std::move(joint)}
    {
    }

    [[nodiscard]] constexpr auto mode_values() noexcept -> mode_topology_type&
    {
        return mode_values_;
    }

    [[nodiscard]] constexpr auto mode_values() const noexcept
        -> const mode_topology_type&
    {
        return mode_values_;
    }

    [[nodiscard]] constexpr auto joint() noexcept -> JointT& { return joint_; }

    [[nodiscard]] constexpr auto joint() const noexcept -> const JointT&
    {
        return joint_;
    }

   private:
    mode_topology_type mode_values_;
    JointT joint_;
};

template <typename ModeTopology, typename JointT>
JointTopology(ModeTopology, JointT) -> JointTopology<ModeTopology, JointT>;

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_JOINT_TOPOLOGY_H_
