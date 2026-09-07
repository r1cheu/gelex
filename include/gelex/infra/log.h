// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_INFRA_LOG_H_
#define GELEX_INFRA_LOG_H_

#include <cstdint>
#include <functional>
#include <string_view>

namespace gelex
{

enum class Level : std::uint8_t
{
    Info,
    Warn,
    Error
};

using Sink = std::function<void(Level, std::string_view)>;

void set_sink(Sink sink);

void info(std::string_view message);
void warn(std::string_view message);
void error(std::string_view message);

}  // namespace gelex

#endif  // GELEX_INFRA_LOG_H_
