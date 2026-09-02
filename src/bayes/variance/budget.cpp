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

#include "gelex/bayes/variance/budget.h"

#include <cmath>
#include <fmt/format.h>
#include <string_view>

#include "gelex/exception.h"

namespace gelex
{

namespace
{

auto validate_variance_share(double share, std::string_view name) -> void
{
    if (!std::isfinite(share) || share < 0.0)
    {
        throw GelexException(
            fmt::format(
                "{} variance share must be finite and non-negative, got {}",
                name,
                share));
    }
}

}  // namespace

VarianceBudget::VarianceBudget(Shares shares)
    : genetic_{shares.additive, shares.dominance}, random_{shares.random}
{
    validate_variance_share(shares.additive, "additive");
    validate_variance_share(shares.dominance, "dominance");
    validate_variance_share(shares.random, "random");

    const double remaining = residual();
    if (!(remaining > 0.0))
    {
        throw GelexException(
            fmt::format(
                "variance shares must sum to less than 1, got {}",
                1.0 - remaining));
    }
}

}  // namespace gelex
