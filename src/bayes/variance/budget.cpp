// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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

VarianceBudget::VarianceBudget(Proportion shares)
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
