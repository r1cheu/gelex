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

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

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

    virtual auto on(
        std::string_view name,
        std::span<const std::string> value,
        FieldFlag flags) -> void = 0;

    virtual auto on(
        std::string_view name,
        std::string_view value,
        FieldFlag flags) -> void = 0;

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

    [[nodiscard]] auto path() const -> std::string_view { return path_; }

    [[nodiscard]] auto field_path(std::string_view name) const -> std::string
    {
        auto result = path_;
        if (!result.empty())
        {
            result.append("/");
        }
        result.append(name);
        return result;
    }

   private:
    auto enter(std::string_view name) -> void
    {
        scope_sizes_.push_back(path_.size());
        if (!path_.empty())
        {
            path_.append("/");
        }
        path_.append(name);
    }

    auto leave() -> void
    {
        path_.resize(scope_sizes_.back());
        scope_sizes_.pop_back();
    }

    std::string path_;
    std::vector<std::size_t> scope_sizes_;
};

}  // namespace gelex::infra

#endif  // GELEX_INFRA_FIELD_VISITOR_H_
