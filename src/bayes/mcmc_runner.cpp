// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/mcmc_runner.h"

#include <fmt/format.h>

#include "gelex/exception.h"

namespace gelex
{

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
MCMCRunner::MCMCRunner(int iterations, int burn_in, int thin)
    : iterations_{iterations}, burn_in_{burn_in}, thin_{thin}
{
    if (iterations_ <= 0)
    {
        throw GelexException(
            fmt::format("iterations must be positive, got {}", iterations_));
    }
    if (burn_in_ < 0 || burn_in_ >= iterations_)
    {
        throw GelexException(
            fmt::format(
                "burn-in must satisfy 0 <= burn-in < iterations, got {} "
                "(iterations={})",
                burn_in_,
                iterations_));
    }
    if (thin_ <= 0)
    {
        throw GelexException(
            fmt::format("thin must be positive, got {}", thin_));
    }
    if ((iterations_ - burn_in_) % thin_ != 0)
    {
        throw GelexException(
            fmt::format(
                "thin ({}) must divide iterations - burn-in ({})",
                thin_,
                iterations_ - burn_in_));
    }
}

}  // namespace gelex
