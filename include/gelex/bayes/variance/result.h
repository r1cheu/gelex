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

#ifndef GELEX_BAYES_VARIANCE_RESULT_H_
#define GELEX_BAYES_VARIANCE_RESULT_H_

#include <utility>

#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/detail/result_factory.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/variance/draws.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

template <GeneticModeSet Modes>
class VarianceSummaryResult
{
   public:
    VarianceSummaryResult(
        HomogeneousModeValues<Modes, ScalarResult> explained_variance,
        HomogeneousModeValues<Modes, ScalarResult> heritability,
        ScalarResult total_explained_variance,
        ScalarResult total_heritability)
        : explained_variance_{std::move(explained_variance)},
          heritability_{std::move(heritability)},
          total_explained_variance_{std::move(total_explained_variance)},
          total_heritability_{std::move(total_heritability)}
    {
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto explained_variance() const noexcept
        -> const ScalarResult&
    {
        return explained_variance_.template get<Mode>();
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto heritability() const noexcept -> const ScalarResult&
    {
        return heritability_.template get<Mode>();
    }

    [[nodiscard]] auto explained_variances() const noexcept
        -> const HomogeneousModeValues<Modes, ScalarResult>&
    {
        return explained_variance_;
    }

    [[nodiscard]] auto heritabilities() const noexcept
        -> const HomogeneousModeValues<Modes, ScalarResult>&
    {
        return heritability_;
    }

    [[nodiscard]] auto total_explained_variance() const noexcept
        -> const ScalarResult&
    {
        return total_explained_variance_;
    }

    [[nodiscard]] auto total_heritability() const noexcept
        -> const ScalarResult&
    {
        return total_heritability_;
    }

   private:
    HomogeneousModeValues<Modes, ScalarResult> explained_variance_;
    HomogeneousModeValues<Modes, ScalarResult> heritability_;
    ScalarResult total_explained_variance_;
    ScalarResult total_heritability_;
};

namespace detail
{

template <GeneticModeSet Modes>
auto make_result(const VarianceSummaryDraws<Modes>& draws)
    -> VarianceSummaryResult<Modes>
{
    auto explained_variance = generate_mode_values<Modes>(
        [&]<GeneticMode Mode>()
        { return make_result(draws.template explained_variance<Mode>()); });
    auto heritability = generate_mode_values<Modes>(
        [&]<GeneticMode Mode>()
        { return make_result(draws.template heritability<Mode>()); });
    auto total_explained_variance
        = make_result(draws.total_explained_variance());
    auto total_heritability = make_result(draws.total_heritability());
    return VarianceSummaryResult<Modes>{
        std::move(explained_variance),
        std::move(heritability),
        std::move(total_explained_variance),
        std::move(total_heritability)};
}

}  // namespace detail

}  // namespace gelex

#endif  // GELEX_BAYES_VARIANCE_RESULT_H_
