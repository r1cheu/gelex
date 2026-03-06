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

#ifndef GELEX_INFRA_DETAIL_INDICATOR_H_
#define GELEX_INFRA_DETAIL_INDICATOR_H_

#include <cstddef>
#include <memory>
#include <string_view>

#include "barkeep.h"

namespace gelex
{
namespace detail
{

inline const barkeep::BarParts BAR_STYLE{
    .left = "[",
    .right = "]",
    .fill = {"\033[1;36m━\033[0m"},
    .empty = {"-"}};

inline const barkeep::Strings GREEN_SPINNER{
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

}  // namespace detail
}  // namespace gelex

#endif  // GELEX_INFRA_DETAIL_INDICATOR_H_
