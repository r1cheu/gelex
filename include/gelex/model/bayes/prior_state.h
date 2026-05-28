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

#ifndef GELEX_MODEL_BAYES_PRIOR_STATE_H_
#define GELEX_MODEL_BAYES_PRIOR_STATE_H_

#include <cstddef>
#include <string_view>

#include <fmt/format.h>
#include <Eigen/Core>
#include <vector>

#include "gelex/exception.h"
#include "gelex/model/bayes/prior_parameters.h"

namespace gelex::infra
{
class FieldVisitor;
}

namespace gelex::bayes
{

struct ProportionState
{
    static constexpr std::string_view name = "proportion";

    ProportionState() = default;
    ProportionState(
        const MixtureProportion& proportion,
        Eigen::Index num_markers);

    auto visit(infra::FieldVisitor& visitor) -> void;

    Eigen::VectorXi assignment;
    Eigen::VectorXi count;
    Eigen::VectorXd value;
    UpdatePolicy update{UpdatePolicy::fixed};
};

struct ComponentState
{
    static constexpr std::string_view name = "component";

    ComponentState() = default;
    ComponentState(Eigen::Index num_components, Eigen::Index num_individuals);

    auto visit(infra::FieldVisitor& visitor) -> void;

    std::vector<Eigen::VectorXd> gebv;
    Eigen::VectorXd gebv_var;
};

class SingleGeneticPriorState
{
   public:
    static constexpr std::string_view name = "single";

    auto operator=(const SingleGeneticPriorState&)
        -> SingleGeneticPriorState& = delete;
    auto operator=(SingleGeneticPriorState&&) noexcept
        -> SingleGeneticPriorState& = delete;

    virtual ~SingleGeneticPriorState() = default;

    template <typename Capability>
    auto get_if() -> Capability*
    {
        return dynamic_cast<Capability*>(this);
    }

    template <typename Capability>
    auto get_if() const -> const Capability*
    {
        return dynamic_cast<const Capability*>(this);
    }

    template <typename Capability>
    auto get() -> Capability&
    {
        auto* capability = get_if<Capability>();
        if (capability == nullptr)
        {
            throw GelexException(
                fmt::format(
                    "single genetic prior state lacks required capability: {}",
                    Capability::name));
        }
        return *capability;
    }

    template <typename Capability>
    auto get() const -> const Capability&
    {
        const auto* capability = get_if<Capability>();
        if (capability == nullptr)
        {
            throw GelexException(
                fmt::format(
                    "single genetic prior state lacks required capability: {}",
                    Capability::name));
        }
        return *capability;
    }

    virtual auto visit(infra::FieldVisitor& visitor) -> void = 0;

   protected:
    SingleGeneticPriorState() = default;
    SingleGeneticPriorState(const SingleGeneticPriorState&) = default;
    SingleGeneticPriorState(SingleGeneticPriorState&&) noexcept = default;
};

class JointGeneticPriorState
{
   public:
    static constexpr std::string_view name = "joint";

    auto operator=(const JointGeneticPriorState&)
        -> JointGeneticPriorState& = delete;
    auto operator=(JointGeneticPriorState&&) noexcept
        -> JointGeneticPriorState& = delete;

    virtual ~JointGeneticPriorState() = default;

    template <typename Capability>
    auto get_if() -> Capability*
    {
        return dynamic_cast<Capability*>(this);
    }

    template <typename Capability>
    auto get_if() const -> const Capability*
    {
        return dynamic_cast<const Capability*>(this);
    }

    template <typename Capability>
    auto get() -> Capability&
    {
        auto* capability = get_if<Capability>();
        if (capability == nullptr)
        {
            throw GelexException(
                fmt::format(
                    "joint genetic prior state lacks required capability: {}",
                    Capability::name));
        }
        return *capability;
    }

    template <typename Capability>
    auto get() const -> const Capability&
    {
        const auto* capability = get_if<Capability>();
        if (capability == nullptr)
        {
            throw GelexException(
                fmt::format(
                    "joint genetic prior state lacks required capability: {}",
                    Capability::name));
        }
        return *capability;
    }

    virtual auto visit(infra::FieldVisitor& visitor) -> void = 0;

   protected:
    JointGeneticPriorState() = default;
    JointGeneticPriorState(const JointGeneticPriorState&) = default;
    JointGeneticPriorState(JointGeneticPriorState&&) noexcept = default;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_PRIOR_STATE_H_
