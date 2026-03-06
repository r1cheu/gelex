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

#include "gelex/infra/detail/indicator.h"

#include <unistd.h>
#include <cstdio>
#include <vector>

#include <barkeep.h>

namespace gelex
{
namespace detail
{
namespace bk = barkeep;

namespace
{

auto should_use_no_tty() -> bool
{
    return isatty(fileno(stdout)) == 0;
}

}  // namespace

auto create_progress_bar(size_t& counter, size_t total, std::string_view format)
    -> ProgressBar
{
    const bool no_tty = should_use_no_tty();

    std::vector<std::shared_ptr<bk::BaseDisplay>> elements;

    elements.push_back((bk::Animation(
        {.message = " ",
         .style = gelex::detail::GREEN_SPINNER,
         .interval = 0.08,
         .no_tty = no_tty,
         .show = false})));
    auto before = bk::Status(
        {.style = bk::Strings{" "}, .no_tty = no_tty, .show = false});
    elements.push_back(before);
    elements.push_back(
        bk::ProgressBar(
            &counter,
            {.total = total,
             .format = std::string(format),
             .speed = 0.1,
             .style = gelex::detail::BAR_STYLE,
             .no_tty = no_tty,
             .show = false}));
    auto after = bk::Status(
        {.style = bk::Strings{" "}, .no_tty = no_tty, .show = false});
    elements.push_back(after);

    return {
        .display = bk::Composite(elements, ""),
        .before_bar = before,
        .after_bar = after};
}

auto create_progress_info() -> ProgressInfo
{
    const bool no_tty = should_use_no_tty();

    std::vector<std::shared_ptr<barkeep::BaseDisplay>> elements;

    elements.push_back(
        barkeep::Animation(
            {.message = " ",
             .style = GREEN_SPINNER,
             .interval = 0.08,
             .no_tty = no_tty,
             .show = false}));
    auto status = barkeep::Status(
        {.message = " ",
         .style = barkeep::Strings{" "},
         .no_tty = no_tty,
         .show = false});
    elements.push_back(status);

    return {
        .display = barkeep::Composite(elements, ""), .progress_info = status};
}

}  // namespace detail
}  // namespace gelex
