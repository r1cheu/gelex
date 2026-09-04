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

#ifndef APPS_CLI_SUMMARY_H_
#define APPS_CLI_SUMMARY_H_

#include <fmt/format.h>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cli/formatter.h"

namespace cli
{

class Summary
{
   public:
    explicit Summary(std::string_view title);

    template <typename... Args>
    auto field(
        std::string_view label,
        fmt::format_string<Args...> fmt_str,
        Args&&... args) -> Summary&
    {
        lines_.push_back(
            cli::field(label, fmt_str, std::forward<Args>(args)...));
        return *this;
    }

    auto show() const -> void;

   private:
    std::vector<std::string> lines_;
};

}  // namespace cli

#endif  // APPS_CLI_SUMMARY_H_
