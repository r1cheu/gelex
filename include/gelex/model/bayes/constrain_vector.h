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

#ifndef GELEX_MODEL_BAYES_CONSTRAIN_VECTOR_H_
#define GELEX_MODEL_BAYES_CONSTRAIN_VECTOR_H_

#include <cmath>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <span>
#include <string_view>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/exception.h"

namespace gelex::bayes
{

template <typename C, typename T>
concept VectorConstraint
    = std::floating_point<T> && requires(std::span<const T> values) {
          { C::name } -> std::convertible_to<std::string_view>;
          { C::template check<T>(values) } -> std::same_as<void>;
      };

namespace detail
{
struct Positive
{
    static constexpr std::string_view name = "PositiveVector";

    template <std::floating_point T>
    static void check(std::span<const T> values)
    {
        for (T v : values)
        {
            if (!(v > T{0}))
            {
                throw GelexException(
                    fmt::format(
                        "{}: all entries must be strictly positive", name));
            }
        }
    }
};

struct Simplex
{
    static constexpr std::string_view name = "Simplex";

    template <std::floating_point T>
    static void check(std::span<const T> values)
    {
        if (values.size() < 2)
        {
            throw GelexException(
                fmt::format(
                    "{}: must have at least 2 components but got {}",
                    name,
                    values.size()));
        }
        T sum{0};
        for (T v : values)
        {
            if (!(v > T{0}))
            {
                throw GelexException(
                    fmt::format(
                        "{}: all entries must be strictly positive", name));
            }
            sum += v;
        }
        const T tol = std::numeric_limits<T>::epsilon()
                      * static_cast<T>(values.size()) * T{64};
        if (std::abs(sum - T{1}) > tol)
        {
            throw GelexException(
                fmt::format("{}: entries must sum to 1", name));
        }
    }
};

struct ScaleMultiplier
{
    static constexpr std::string_view name = "ScaleMultiplier";

    template <std::floating_point T>
    static void check(std::span<const T> values)
    {
        if (values.size() < 2)
        {
            throw GelexException(
                fmt::format("{}: must have at least 2 elements", name));
        }
        if (values[0] != T{0})
        {
            throw GelexException(
                fmt::format("{}: first element must be 0", name));
        }
        for (std::size_t i = 1; i < values.size(); ++i)
        {
            if (!std::isfinite(values[i]))
            {
                throw GelexException(
                    fmt::format("{}: element {} must be finite", name, i));
            }
            if (!(values[i] > T{0}))
            {
                throw GelexException(
                    fmt::format("{}: element {} must be positive", name, i));
            }
        }
    }
};

}  // namespace detail

template <std::floating_point T, typename Constraint>
    requires VectorConstraint<Constraint, T>
class ConstrainVector
{
   public:
    ConstrainVector(std::initializer_list<T> values)
        : ConstrainVector(std::span<const T>{values.begin(), values.size()})
    {
    }

    explicit ConstrainVector(std::span<const T> values)
        : values_(values.begin(), values.end())
    {
        if (values_.empty())
        {
            throw GelexException(
                fmt::format(
                    "{}: must contain at least one element", Constraint::name));
        }
        Constraint::template check<T>(values_);
    }

    [[nodiscard]] auto size() const noexcept -> std::size_t
    {
        return values_.size();
    }

    [[nodiscard]] auto data() const noexcept -> const T*
    {
        return values_.data();
    }

    [[nodiscard]] auto operator[](std::size_t i) const -> T
    {
        return values_[i];
    }

    [[nodiscard]] auto begin() const noexcept { return values_.begin(); }
    [[nodiscard]] auto end() const noexcept { return values_.end(); }

    [[nodiscard]] auto to_vector() const -> const std::vector<T>&
    {
        return values_;
    }

    [[nodiscard]] auto to_mat() const -> Eigen::VectorX<T>
    {
        return Eigen::Map<const Eigen::VectorX<T>>(
            values_.data(), static_cast<Eigen::Index>(values_.size()));
    }

   private:
    std::vector<T> values_;
};

template <std::floating_point T>
using Simplex = ConstrainVector<T, detail::Simplex>;

template <std::floating_point T>
using PositiveVector = ConstrainVector<T, detail::Positive>;

template <std::floating_point T>
using ScaleMultiplier = ConstrainVector<T, detail::ScaleMultiplier>;

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_CONSTRAIN_VECTOR_H_
