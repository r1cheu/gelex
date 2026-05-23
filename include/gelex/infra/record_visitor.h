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

#ifndef GELEX_INFRA_RECORD_VISITOR_H_
#define GELEX_INFRA_RECORD_VISITOR_H_

#include <cstddef>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include <fmt/format.h>
#include <Eigen/Core>

namespace gelex::infra
{

template <bool Mutable>
class BasicRecordVisitor
{
   public:
    template <typename Value>
    using RecordType = std::conditional_t<
        Mutable,
        Eigen::Ref<Value>,
        const Eigen::Ref<const Value>&>;

    virtual ~BasicRecordVisitor() = default;

    virtual auto visit(std::string_view path, RecordType<Eigen::VectorXf> value)
        -> void
        = 0;
    virtual auto visit(std::string_view path, RecordType<Eigen::VectorXd> value)
        -> void
        = 0;
    virtual auto visit(std::string_view path, RecordType<Eigen::VectorXi> value)
        -> void
        = 0;

    template <typename Value>
    auto emit(
        std::string_view capability,
        std::size_t slot,
        std::string_view field,
        Value&& value) -> void
    {
        visit(make_path(capability, slot, field), std::forward<Value>(value));
    }

   protected:
    BasicRecordVisitor() = default;
    BasicRecordVisitor(const BasicRecordVisitor&) = default;
    auto operator=(const BasicRecordVisitor&) -> BasicRecordVisitor& = default;
    BasicRecordVisitor(BasicRecordVisitor&&) noexcept = default;
    auto operator=(BasicRecordVisitor&&) noexcept
        -> BasicRecordVisitor& = default;

   private:
    static auto make_path(
        std::string_view capability,
        std::size_t slot,
        std::string_view field) -> std::string
    {
        return fmt::format("{}/{}/{}", capability, slot, field);
    }
};

using RecordVisitor = BasicRecordVisitor<false>;
using MutableRecordVisitor = BasicRecordVisitor<true>;

}  // namespace gelex::infra

#endif  // GELEX_INFRA_RECORD_VISITOR_H_
