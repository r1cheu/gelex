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

#ifndef APPS_CLI_THEME_H_
#define APPS_CLI_THEME_H_

#include <cstdint>
#include <fmt/color.h>
#include <string>

namespace cli
{

// Semantic color roles. Each maps to one of the 16 ANSI named colors, which the
// terminal renders using the user's configured theme palette rather than a
// fixed RGB value, so output follows whatever theme the shell is set to.
enum class ColorRole : std::uint8_t
{
    heading,
    option_name,
    option_value,
    success,
    warning,
    error,
    accent,
    muted,
};

// Whether ANSI styling should be emitted. Honors the NO_COLOR and
// CLICOLOR_FORCE conventions on top of TTY detection: NO_COLOR disables,
// CLICOLOR_FORCE forces, otherwise falls back to is_tty().
auto should_colorize() -> bool;

auto style_for(ColorRole role) -> fmt::text_style;

auto colorize(ColorRole role, std::string text) -> std::string;

auto colorize(const fmt::text_style& style, std::string text) -> std::string;

// "[error] " prefix with a themed red "error"; gated by should_colorize().
auto error_marker() -> std::string;

}  // namespace cli

#endif  // APPS_CLI_THEME_H_
