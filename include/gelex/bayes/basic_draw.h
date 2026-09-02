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

#ifndef GELEX_BAYES_BASIC_DRAW_H_
#define GELEX_BAYES_BASIC_DRAW_H_

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <limits>
#include <span>
#include <string_view>
#include <utility>

#include "gelex/bayes/stats/detail/running_stats.h"
#include "gelex/bayes/stats/result.h"
#include "gelex/exception.h"
#include "gelex/io/binary_writer.h"

namespace gelex
{

namespace detail
{

inline auto to_eigen_draw_size(std::uint64_t rows) -> Eigen::Index
{
    if (!std::in_range<Eigen::Index>(rows))
    {
        throw GelexException("draw row count is out of range");
    }
    return static_cast<Eigen::Index>(rows);
}

inline auto throw_if_empty(const auto& stats, std::string_view identifier)
    -> void
{
    if (stats.empty())
    {
        throw GelexException(
            fmt::format("posterior '{}' has no recorded draws", identifier));
    }
}

}  // namespace detail

struct EmptyDraw
{
    auto append(const auto& /*ignored*/) const noexcept -> void {}
};

class ScalarDraw
{
   public:
    explicit ScalarDraw(PayloadWriter<double> payload) noexcept
        : payload_(std::move(payload))
    {
    }

    auto append(double value) -> void
    {
        payload_.append(value);
        stats_.update(value);
    }

    [[nodiscard]] auto result() const -> ScalarRunningStatsResult
    {
        detail::throw_if_empty(stats_, identifier());
        return stats_.result();
    }

    [[nodiscard]] auto identifier() const noexcept -> std::string_view
    {
        return payload_.identifier();
    }

   private:
    PayloadWriter<double> payload_;
    detail::ScalarRunningStats stats_;
};

// Traces are stored as float32: their precision is bounded by Monte Carlo
// error, while the running statistics stay in double.
class VectorDraw
{
   public:
    explicit VectorDraw(PayloadWriter<float> payload)
        : stats_(detail::to_eigen_draw_size(payload.rows())),
          buffer_(detail::to_eigen_draw_size(payload.rows())),
          payload_(std::move(payload))
    {
    }

    auto append(std::span<const double> values) -> void
    {
        const Eigen::Map<const Eigen::VectorXd> view{
            values.data(), static_cast<Eigen::Index>(values.size())};
        stats_.update(view);
        buffer_ = view.cast<float>();
        payload_.append(buffer_);
    }

    [[nodiscard]] auto result() const -> VectorRunningStatsResult
    {
        detail::throw_if_empty(stats_, identifier());
        return stats_.result();
    }

    [[nodiscard]] auto mean_square() const -> Eigen::VectorXd
    {
        detail::throw_if_empty(stats_, identifier());
        return stats_.mean_square();
    }

    [[nodiscard]] auto identifier() const noexcept -> std::string_view
    {
        return payload_.identifier();
    }

   private:
    detail::VectorRunningStats stats_;
    Eigen::VectorXf buffer_;
    PayloadWriter<float> payload_;
};

template <std::size_t CategoryCount>
    requires(
        CategoryCount > 1
        && CategoryCount <= static_cast<std::size_t>(
                                std::numeric_limits<std::uint8_t>::max())
                                + 1)
class CategoryDraw
{
   public:
    explicit CategoryDraw(PayloadWriter<std::uint8_t> payload)
        : stats_(detail::to_eigen_draw_size(payload.rows())),
          payload_(std::move(payload))
    {
    }

    auto append(std::span<const std::uint8_t> values) -> void
    {
        payload_.append(values);
        stats_.update(values);
    }

    [[nodiscard]] auto result() const -> CategoryRunningStatsResult
    {
        detail::throw_if_empty(stats_, identifier());
        return stats_.result();
    }

    template <typename IsIncluded>
    [[nodiscard]] auto probability_of(IsIncluded is_included) const
        -> Eigen::VectorXd
    {
        detail::throw_if_empty(stats_, identifier());
        return stats_.probability_of(std::move(is_included));
    }

    [[nodiscard]] auto identifier() const noexcept -> std::string_view
    {
        return payload_.identifier();
    }

   private:
    detail::CategoryRunningStats<CategoryCount> stats_;
    PayloadWriter<std::uint8_t> payload_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_BASIC_DRAW_H_
