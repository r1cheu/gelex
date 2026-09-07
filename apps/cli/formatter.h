// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_FORMATTER_H_
#define APPS_CLI_FORMATTER_H_

#include <Eigen/Core>
#include <cstddef>
#include <fmt/base.h>
#include <fmt/color.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <string>
#include <string_view>
#include <utility>

#include "theme.h"

namespace cli
{

std::string command_banner(
    std::string_view version,
    std::string_view task,
    size_t width = 80);

// If s[pos..] begins an ANSI CSI escape (ESC '[' ... final byte 0x40-0x7e),
// returns the index just past it; otherwise returns pos unchanged. pos must be
// a valid index into s.
size_t skip_ansi_escape(std::string_view s, size_t pos) noexcept;

std::string format_eta(double seconds);
std::string done_message(double elapsed_seconds);

template <typename... Args>
std::string section(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return " "
           + colorize(
               style_for(ColorRole::accent) | fmt::emphasis::bold,
               fmt::format(fmt_str, std::forward<Args>(args)...));
}

template <typename... Args>
std::string success(fmt::format_string<Args...> fmt_str, Args&&... args)
{
    return colorize(ColorRole::success, "  ✓ ")
           + fmt::format(fmt_str, std::forward<Args>(args)...);
}

// Unified completion line shared by every command: "✓ Results saved:
// <prefix>*  (<artifacts>)". prefix is the output path stem (the trailing '*'
// makes it directly glob-able across shells); artifacts lists the produced
// suffixes, always including .log.
inline std::string results_saved(
    std::string_view prefix,
    std::string_view artifacts)
{
    return success("Results saved: {}*  ({})", prefix, artifacts);
}

inline constexpr size_t field_label_width = 16;

// Aligned "   <label>: <value>" line for summary sections; single source of the
// label column width so the alignment stays consistent across reporters.
template <typename... Args>
std::string field(
    std::string_view label,
    fmt::format_string<Args...> fmt_str,
    Args&&... args)
{
    return fmt::format(
        "   {:<{}}: {}",
        label,
        field_label_width,
        fmt::format(fmt_str, std::forward<Args>(args)...));
}

struct AbbrNumber
{
    explicit AbbrNumber(size_t value) : value(static_cast<double>(value)) {}
    explicit AbbrNumber(Eigen::Index value) : value(static_cast<double>(value))
    {
    }
    explicit AbbrNumber(double value) : value(value) {}
    double value;
};

}  // namespace cli

namespace fmt
{
template <>
struct formatter<cli::AbbrNumber> : formatter<double>
{
    auto format(cli::AbbrNumber hr, format_context& ctx) const
        -> format_context::iterator;
};
}  // namespace fmt

#endif  // APPS_CLI_FORMATTER_H_
