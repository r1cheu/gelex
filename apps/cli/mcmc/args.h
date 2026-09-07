// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_MCMC_ARGS_H_
#define APPS_CLI_MCMC_ARGS_H_

namespace CLI
{
class App;
}

auto setup_mcmc_command(CLI::App& program, int& exit_code) -> void;

#endif  // APPS_CLI_MCMC_ARGS_H_
