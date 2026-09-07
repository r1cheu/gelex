// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
