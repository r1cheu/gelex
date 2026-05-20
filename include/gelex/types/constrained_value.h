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

#ifndef GELEX_TYPES_CONSTRAINED_VALUE_H_
#define GELEX_TYPES_CONSTRAINED_VALUE_H_

#include <cmath>
#include <concepts>
#include <string_view>

#include <fmt/format.h>

#include "gelex/exception.h"

namespace gelex
{

template <typename C, typename T>
concept ValueConstraint = std::floating_point<T> && requires(T value) {
    { C::name } -> std::convertible_to<std::string_view>;
    { C::template check<T>(value) } -> std::same_as<void>;
};

namespace detail
{

struct PositiveScalar
{
    static constexpr std::string_view name = "PositiveScalar";

    template <std::floating_point T>
    static void check(T value)
    {
        if (!std::isfinite(value) || !(value > T{0}))
        {
            throw GelexException(
                fmt::format(
                    "{}: value must be finite and strictly positive but got {}",
                    name,
                    value));
        }
    }
};

struct OpenUnitInterval
{
    static constexpr std::string_view name = "OpenUnitInterval";

    template <std::floating_point T>
    static void check(T value)
    {
        if (!std::isfinite(value) || !(value > T{0}) || !(value < T{1}))
        {
            throw GelexException(
                fmt::format(
                    "{}: value must be in (0, 1) but got {}", name, value));
        }
    }
};

}  // namespace detail

template <std::floating_point T, typename Constraint>
    requires ValueConstraint<Constraint, T>
class ConstrainedValue
{
   public:
    explicit ConstrainedValue(T value) : value_(value)
    {
        Constraint::template check<T>(value_);
    }

    [[nodiscard]] auto value() const noexcept -> T { return value_; }

    explicit operator T() const { return value_; }

   private:
    T value_;
};

template <std::floating_point T>
using PositiveScalar = ConstrainedValue<T, detail::PositiveScalar>;

template <std::floating_point T>
using OpenUnitInterval = ConstrainedValue<T, detail::OpenUnitInterval>;

}  // namespace gelex

#endif  // GELEX_TYPES_CONSTRAINED_VALUE_H_
