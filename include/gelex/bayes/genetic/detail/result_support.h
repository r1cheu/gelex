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

#ifndef GELEX_BAYES_GENETIC_DETAIL_RESULT_SUPPORT_H_
#define GELEX_BAYES_GENETIC_DETAIL_RESULT_SUPPORT_H_

#include <span>
#include <string>
#include <vector>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/genetic/traits.h"
#include "gelex/bayes/genetic_family.h"

namespace gelex::detail
{

inline auto make_result(const EmptyDraw& /*draw*/) -> EmptyPosteriorResult
{
    return {};
}

inline auto make_result(const ScalarDraw& draw) -> ScalarPosteriorResult
{
    return ScalarPosteriorResult{std::string{draw.identifier()}, draw.result()};
}

inline auto make_result(const VectorDraw& draw) -> VectorPosteriorResult
{
    return VectorPosteriorResult{std::string{draw.identifier()}, draw.result()};
}

inline auto make_result(
    const VectorDraw& draw,
    std::span<const std::string> column_names) -> CoefficientPosteriorResult
{
    return CoefficientPosteriorResult{
        make_result(draw),
        std::vector<std::string>{column_names.begin(), column_names.end()}};
}

template <VarianceLayout Kind>
auto make_marker_variance_result(const marker_variance_draw_t<Kind>& draw)
    -> marker_variance_result_t<Kind>
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return make_result(draw);
    }
    else
    {
        static_cast<void>(draw);
        return {};
    }
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_RESULT_SUPPORT_H_
