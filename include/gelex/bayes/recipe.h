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

#ifndef GELEX_BAYES_RECIPE_H_
#define GELEX_BAYES_RECIPE_H_

#include <fmt/format.h>
#include <utility>

#include "gelex/bayes/variance/budget.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

template <GeneticModeSet Modes, typename GeneticSpec>
class BayesRecipe
{
   public:
    static constexpr GeneticModeSet modes = Modes;
    using genetic_spec_type = GeneticSpec;

    BayesRecipe(genetic_spec_type genetic_spec, VarianceBudget variance)
        : genetic_spec_(std::move(genetic_spec)), variance_(variance)
    {
        validate();
    }

    explicit BayesRecipe(VarianceBudget variance)
        : BayesRecipe(genetic_spec_type{}, variance)
    {
    }

    static auto defaults() -> BayesRecipe
    {
        return BayesRecipe{VarianceBudget{default_proportion(Modes)}};
    }

    auto genetic_spec() const noexcept -> const genetic_spec_type&
    {
        return genetic_spec_;
    }

    auto variance() const noexcept -> const VarianceBudget&
    {
        return variance_;
    }

   private:
    auto validate() const -> void
    {
        if constexpr (requires { GeneticSpec::modes; })
        {
            static_assert(
                GeneticSpec::modes == Modes,
                "genetic spec modes must match recipe modes");
        }

        for (const auto mode : all_genetic_modes)
        {
            const double proportion = variance_.genetic(mode);
            const bool is_present = Modes.contains(mode);
            if (is_present && proportion == 0.0)
            {
                throw GelexException(
                    fmt::format(
                        "invalid Bayes recipe input: {} variance proportion "
                        "must be positive when the mode is present, got {}",
                        mode,
                        proportion));
            }
            if (!is_present && proportion != 0.0)
            {
                throw GelexException(
                    fmt::format(
                        "invalid Bayes recipe input: {} variance proportion "
                        "must be zero when the mode is absent, got {}",
                        mode,
                        proportion));
            }
        }
    }

    [[no_unique_address]] genetic_spec_type genetic_spec_;
    VarianceBudget variance_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_RECIPE_H_
