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

#include "gelex/bayes/detail/recipe_validation.h"

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <string>
#include <vector>

#include "gelex/bayes/variance_budget.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

namespace
{

class RecipeIssues
{
   public:
    auto add(GeneticModeSet scope, std::string issue) -> void
    {
        issues_.push_back(fmt::format("{} {}", scope, issue));
    }

    auto throw_if_any() const -> void
    {
        if (issues_.empty())
        {
            return;
        }

        throw GelexException(
            fmt::format(
                "invalid Bayes recipe input:\n  {}",
                fmt::join(issues_, "\n  ")));
    }

   private:
    std::vector<std::string> issues_;
};

auto check(
    RecipeIssues& issues,
    const VarianceBudget& budget,
    GeneticModeSet modes) -> void
{
    for (const auto mode : all_genetic_modes)
    {
        const auto share = budget.share(mode);
        const auto scope = GeneticModeSet{mode};

        if (modes.contains(mode))
        {
            if (share == 0.0)
            {
                issues.add(
                    scope,
                    fmt::format(
                        "variance share must be positive when the mode is "
                        "present, got {}",
                        share));
            }
        }
        else if (share != 0.0)
        {
            issues.add(
                scope,
                fmt::format(
                    "variance share must be zero when the mode is absent, "
                    "got {}",
                    share));
        }
    }
}

}  // namespace

auto validate_recipe_inputs(
    const VarianceBudget& variance,
    GeneticModeSet modes) -> void
{
    auto issues = RecipeIssues{};
    check(issues, variance, modes);
    issues.throw_if_any();
}

}  // namespace gelex::detail
