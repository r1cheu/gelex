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

#ifndef GELEX_BAYES_BASIC_RESULT_H_
#define GELEX_BAYES_BASIC_RESULT_H_

#include <fmt/format.h>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/stats/result.h"
#include "gelex/exception.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

struct EmptyResult
{
};

class ScalarResult
{
   public:
    ScalarResult(std::string identifier, ScalarRunningStatsResult statistics)
        : identifier_{std::move(identifier)}, statistics_{statistics}
    {
        if (identifier_.empty())
        {
            throw GelexException("posterior identifier is empty");
        }
    }

    [[nodiscard]] auto identifier() const noexcept -> std::string_view
    {
        return identifier_;
    }

    [[nodiscard]] auto statistics() const noexcept
        -> const ScalarRunningStatsResult&
    {
        return statistics_;
    }

   private:
    std::string identifier_;
    ScalarRunningStatsResult statistics_;
};

class VectorResult
{
   public:
    VectorResult(std::string identifier, VectorRunningStatsResult statistics)
        : identifier_{std::move(identifier)}, statistics_{std::move(statistics)}
    {
        if (identifier_.empty())
        {
            throw GelexException("posterior identifier is empty");
        }
        if (statistics_.mean.size() != statistics_.stddev.size())
        {
            throw GelexException(
                fmt::format(
                    "posterior '{}' has {} means but {} standard deviations",
                    identifier_,
                    statistics_.mean.size(),
                    statistics_.stddev.size()));
        }
    }

    [[nodiscard]] auto identifier() const noexcept -> std::string_view
    {
        return identifier_;
    }

    [[nodiscard]] auto statistics() const noexcept
        -> const VectorRunningStatsResult&
    {
        return statistics_;
    }

   private:
    std::string identifier_;
    VectorRunningStatsResult statistics_;
};

class CoefficientResult
{
   public:
    CoefficientResult(
        VectorResult posterior,
        std::vector<std::string> column_names)
        : posterior_{std::move(posterior)},
          column_names_{std::move(column_names)}
    {
        if (!std::cmp_equal(
                column_names_.size(), posterior_.statistics().mean.size()))
        {
            throw GelexException(
                fmt::format(
                    "posterior '{}' has {} statistics but {} column names",
                    posterior_.identifier(),
                    posterior_.statistics().mean.size(),
                    column_names_.size()));
        }
    }

    [[nodiscard]] auto identifier() const noexcept -> std::string_view
    {
        return posterior_.identifier();
    }

    [[nodiscard]] auto column_names() const noexcept
        -> std::span<const std::string>
    {
        return column_names_;
    }

    [[nodiscard]] auto statistics() const noexcept
        -> const VectorRunningStatsResult&
    {
        return posterior_.statistics();
    }

   private:
    VectorResult posterior_;
    std::vector<std::string> column_names_;
};

GELEX_NAMESPACE_BEGIN(detail)
inline auto make_result(const EmptyDraw& /*draw*/) -> EmptyResult
{
    return {};
}

inline auto make_result(const ScalarDraw& draw) -> ScalarResult
{
    return ScalarResult{std::string{draw.identifier()}, draw.result()};
}

inline auto make_result(const VectorDraw& draw) -> VectorResult
{
    return VectorResult{std::string{draw.identifier()}, draw.result()};
}

inline auto make_result(
    const VectorDraw& draw,
    std::span<const std::string> column_names) -> CoefficientResult
{
    return CoefficientResult{
        make_result(draw),
        std::vector<std::string>{column_names.begin(), column_names.end()}};
}
GELEX_NAMESPACE_END(detail)

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_BASIC_RESULT_H_
