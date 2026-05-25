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

#ifndef GELEX_INFRA_FIELD_VISITOR_H_
#define GELEX_INFRA_FIELD_VISITOR_H_

#include <string_view>
#include <type_traits>
#include <utility>

#include <Eigen/Core>

#include "gelex/infra/field_flag.h"

namespace gelex::infra
{

class FieldVisitor
{
   public:
    class Scope
    {
       public:
        Scope(FieldVisitor& visitor, std::string_view name) : visitor_(&visitor)
        {
            visitor_->enter(name);
        }

        Scope(const Scope&) = delete;
        Scope(Scope&& other) noexcept
            : visitor_(std::exchange(other.visitor_, nullptr))
        {
        }

        auto operator=(const Scope&) -> Scope& = delete;
        auto operator=(Scope&&) -> Scope& = delete;

        ~Scope()
        {
            if (visitor_ != nullptr)
            {
                visitor_->leave();
            }
        }

       private:
        FieldVisitor* visitor_;
    };

    FieldVisitor(const FieldVisitor&) = delete;
    FieldVisitor(FieldVisitor&&) = delete;
    auto operator=(const FieldVisitor&) -> FieldVisitor& = delete;
    auto operator=(FieldVisitor&&) -> FieldVisitor& = delete;
    virtual ~FieldVisitor() = default;

    [[nodiscard]] auto scope(std::string_view name) -> Scope
    {
        return Scope{*this, name};
    }

    virtual auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXf> value,
        FieldFlag flags) -> void = 0;

    virtual auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXd> value,
        FieldFlag flags) -> void = 0;

    virtual auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXi> value,
        FieldFlag flags) -> void = 0;

    virtual auto on(std::string_view name, double& value, FieldFlag flags)
        -> void = 0;

    virtual auto on(std::string_view name, int& value, FieldFlag flags) -> void
        = 0;

    template <typename Enum>
        requires std::is_enum_v<Enum>
    auto on(std::string_view name, Enum& value, FieldFlag flags) -> void
    {
        auto encoded_value = static_cast<int>(value);
        on(name, encoded_value, flags);
        value = static_cast<Enum>(encoded_value);
    }

   protected:
    FieldVisitor() = default;

   private:
    virtual auto enter(std::string_view name) -> void = 0;
    virtual auto leave() -> void = 0;
};

}  // namespace gelex::infra

#endif  // GELEX_INFRA_FIELD_VISITOR_H_
