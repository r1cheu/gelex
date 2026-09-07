// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_LOGGING_H_
#define APPS_CLI_LOGGING_H_

#include <memory>
#include <spdlog/logger.h>
#include <string_view>

namespace cli::logging
{
void initialize(std::string_view output_prefix = "output");

std::shared_ptr<spdlog::logger>& get();
}  // namespace cli::logging

#endif  // APPS_CLI_LOGGING_H_
