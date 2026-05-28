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

#ifndef GELEX_MODEL_BAYES_GENETIC_PRIOR_H_
#define GELEX_MODEL_BAYES_GENETIC_PRIOR_H_

#include <memory>
#include <string_view>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

class SingleGeneticPrior
{
   public:
    static constexpr std::string_view name = "single";

    auto operator=(const SingleGeneticPrior&) -> SingleGeneticPrior& = delete;
    auto operator=(SingleGeneticPrior&&) noexcept
        -> SingleGeneticPrior& = delete;

    virtual ~SingleGeneticPrior() = default;

    virtual auto mode() const -> GeneticMode = 0;
    virtual auto visit(infra::FieldVisitor& visitor) -> void = 0;

    virtual auto make_state(
        Eigen::Index num_markers,
        Eigen::Index num_individuals) const
        -> std::unique_ptr<SingleGeneticPriorState> = 0;

    template <typename Capability>
    auto get_if() const -> const Capability*
    {
        return dynamic_cast<const Capability*>(this);
    }

    template <typename Capability>
    auto get() const -> const Capability&
    {
        const auto* capability = get_if<Capability>();
        if (capability == nullptr)
        {
            throw GelexException(
                fmt::format(
                    "single genetic prior lacks required capability: {}",
                    Capability::name));
        }
        return *capability;
    }

   protected:
    SingleGeneticPrior() = default;
    SingleGeneticPrior(const SingleGeneticPrior&) = default;
    SingleGeneticPrior(SingleGeneticPrior&&) noexcept = default;
};

class JointGeneticPrior
{
   public:
    static constexpr std::string_view name = "joint";

    auto operator=(const JointGeneticPrior&) -> JointGeneticPrior& = delete;
    auto operator=(JointGeneticPrior&&) noexcept -> JointGeneticPrior& = delete;

    virtual ~JointGeneticPrior() = default;

    virtual auto visit(infra::FieldVisitor& visitor) -> void = 0;

    virtual auto make_state(
        Eigen::Index num_markers,
        Eigen::Index num_individuals) const
        -> std::unique_ptr<JointGeneticPriorState> = 0;

    template <typename Capability>
    auto get_if() const -> const Capability*
    {
        return dynamic_cast<const Capability*>(this);
    }

    template <typename Capability>
    auto get() const -> const Capability&
    {
        const auto* capability = get_if<Capability>();
        if (capability == nullptr)
        {
            throw GelexException(
                fmt::format(
                    "joint genetic prior lacks required capability: {}",
                    Capability::name));
        }
        return *capability;
    }

   protected:
    JointGeneticPrior() = default;
    JointGeneticPrior(const JointGeneticPrior&) = default;
    JointGeneticPrior(JointGeneticPrior&&) noexcept = default;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GENETIC_PRIOR_H_
