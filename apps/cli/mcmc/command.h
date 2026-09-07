// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_MCMC_COMMAND_H_
#define APPS_CLI_MCMC_COMMAND_H_

#include "config.h"

auto mcmc_execute(const cli::McmcConfig& config) -> int;

#endif  // APPS_CLI_MCMC_COMMAND_H_
