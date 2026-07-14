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

#include "theme.h"

#include <cstdlib>
#include <fmt/color.h>
#include <fmt/format.h>
#include <string>
#include <string_view>
#include <utility>

#include "cli/runtime.h"

namespace cli
{

namespace
{

auto env_present(const char* name) -> bool
{
    const char* value = std::getenv(name);
    return value != nullptr && value[0] != '\0';
}

}  // namespace

auto should_colorize() -> bool
{
    if (env_present("NO_COLOR"))
    {
        return false;
    }
    if (const char* force = std::getenv("CLICOLOR_FORCE");
        force != nullptr && force[0] != '\0' && std::string_view{force} != "0")
    {
        return true;
    }
    return is_tty();
}

auto style_for(ColorRole role) -> fmt::text_style
{
    switch (role)
    {
        case ColorRole::heading:
            return fmt::emphasis::bold | fmt::emphasis::underline
                   | fmt::fg(fmt::terminal_color::green);
        case ColorRole::option_name:
            return fmt::fg(fmt::terminal_color::cyan);
        case ColorRole::option_value:
            return fmt::fg(fmt::terminal_color::yellow);
        case ColorRole::success:
            return fmt::fg(fmt::terminal_color::green);
        case ColorRole::warning:
            return fmt::fg(fmt::terminal_color::yellow);
        case ColorRole::error:
            return fmt::fg(fmt::terminal_color::red);
        case ColorRole::accent:
            return fmt::fg(fmt::terminal_color::cyan);
        case ColorRole::muted:
            return fmt::fg(fmt::terminal_color::bright_black);
    }
    return {};
}

auto colorize(const fmt::text_style& style, std::string text) -> std::string
{
    if (text.empty() || !should_colorize())
    {
        return text;
    }
    return fmt::format(style, "{}", text);
}

auto colorize(ColorRole role, std::string text) -> std::string
{
    return colorize(style_for(role), std::move(text));
}

auto error_marker() -> std::string
{
    return fmt::format("[{}] ", colorize(ColorRole::error, "error"));
}

}  // namespace cli
