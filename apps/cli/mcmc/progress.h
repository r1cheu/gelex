// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_MCMC_PROGRESS_H_
#define APPS_CLI_MCMC_PROGRESS_H_

#include <cstddef>

#include "cli/progress.h"

namespace cli
{

class GenotypeProgress
{
   public:
    explicit GenotypeProgress(std::size_t total);
    auto operator()(std::size_t current) -> void;
    auto finish() -> void;

   private:
    cli::Progress progress_;
    decltype(cli::make_rate()) estimate_rate_;
    decltype(cli::make_eta(std::size_t{})) estimate_eta_;
};

class McmcProgress
{
   public:
    McmcProgress(std::size_t total, int burn_in);
    auto operator()(std::size_t current) -> void;
    auto finish() -> void;

   private:
    std::size_t burn_in_;
    cli::Progress progress_;
    decltype(cli::make_rate()) estimate_rate_;
    decltype(cli::make_eta(std::size_t{})) estimate_eta_;
};

}  // namespace cli

#endif  // APPS_CLI_MCMC_PROGRESS_H_
