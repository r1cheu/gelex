// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_RUNTIME_H_
#define APPS_CLI_RUNTIME_H_

namespace cli
{
auto is_tty() -> bool;

auto setup_parallelization(int num_threads) -> void;
}  // namespace cli

#endif  // APPS_CLI_RUNTIME_H_
