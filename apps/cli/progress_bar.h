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

#ifndef APPS_CLI_PROGRESS_BAR_H_
#define APPS_CLI_PROGRESS_BAR_H_

#include <cstddef>
#include <memory>
#include <string_view>

#include "barkeep.h"

namespace cli
{

inline const barkeep::Strings green_spinner{
    "\033[32m⠁\033[0m", "\033[32m⠁\033[0m", "\033[32m⠉\033[0m",
    "\033[32m⠙\033[0m", "\033[32m⠚\033[0m", "\033[32m⠒\033[0m",
    "\033[32m⠂\033[0m", "\033[32m⠂\033[0m", "\033[32m⠒\033[0m",
    "\033[32m⠲\033[0m", "\033[32m⠴\033[0m", "\033[32m⠤\033[0m",
    "\033[32m⠄\033[0m", "\033[32m⠄\033[0m", "\033[32m⠤\033[0m",
    "\033[32m⠠\033[0m", "\033[32m⠠\033[0m", "\033[32m⠤\033[0m",
    "\033[32m⠦\033[0m", "\033[32m⠖\033[0m", "\033[32m⠒\033[0m",
    "\033[32m⠐\033[0m", "\033[32m⠐\033[0m", "\033[32m⠒\033[0m",
    "\033[32m⠓\033[0m", "\033[32m⠋\033[0m", "\033[32m⠉\033[0m",
    "\033[32m⠈\033[0m", "\033[32m⠈\033[0m", "\033[32m \033[0m"};

struct ProgressBar
{
    std::shared_ptr<barkeep::CompositeDisplay> display;
    std::shared_ptr<barkeep::StatusDisplay> before_bar;
    std::shared_ptr<barkeep::StatusDisplay> after_bar;
};

struct ProgressInfo
{
    std::shared_ptr<barkeep::CompositeDisplay> display;
    std::shared_ptr<barkeep::StatusDisplay> progress_info;
};

auto create_progress_bar(
    size_t& counter,
    size_t total,
    std::string_view format = "{bar}") -> ProgressBar;

auto create_progress_info() -> ProgressInfo;

// Erases the just-finished bar line (the one barkeep left with a trailing
// newline) so a persistent summary can take its place. No-op when output is not
// a tty.
auto clear_finished_line() -> void;

}  // namespace cli

#endif  // APPS_CLI_PROGRESS_BAR_H_
