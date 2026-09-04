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

#include "progress.h"

#include <array>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fmt/color.h>
#include <fmt/format.h>
#include <iterator>
#include <mutex>
#include <string>
#include <string_view>
#include <thread>
#include <utility>

#include "cli/formatter.h"
#include "cli/runtime.h"
#include "cli/theme.h"

namespace cli
{
namespace
{

constexpr std::array<std::string_view, 30> spinner_frames{
    "⠁", "⠁", "⠉", "⠙", "⠚", "⠒", "⠂", "⠂", "⠒", "⠲", "⠴", "⠤", "⠄", "⠄", "⠤",
    "⠠", "⠠", "⠤", "⠦", "⠖", "⠒", "⠐", "⠐", "⠒", "⠓", "⠋", "⠉", "⠈", "⠈", " "};
constexpr auto refresh_interval = std::chrono::milliseconds{80};

auto append_bar(
    std::string& output,
    std::size_t current,
    std::size_t total,
    std::size_t width = 20) -> void
{
    std::size_t filled_width = 0;
    if (total != 0)
    {
        filled_width = current >= total ? width : current * width / total;
    }

    if (filled_width != 0)
    {
        if (should_colorize())
        {
            fmt::format_to(
                std::back_inserter(output),
                style_for(ColorRole::accent) | fmt::emphasis::bold,
                "{:━<{}}",
                "",
                filled_width);
        }
        else
        {
            fmt::format_to(
                std::back_inserter(output), "{:━<{}}", "", filled_width);
        }
    }
    output.append(width - filled_width, '-');
}

}  // namespace

Progress::Progress(std::string prefix, std::size_t total, std::string rate_unit)
    : prefix_{std::move(prefix)},
      rate_unit_{std::move(rate_unit)},
      total_{total}
{
    if (std::getenv("GELEX_NO_TTY") == nullptr && is_tty())
    {
        renderer_ = std::jthread{[this](const std::stop_token& stop_token)
                                 { render(stop_token); }};
    }
}

Progress::~Progress()
{
    finish();
}

auto Progress::update(const ProgressSnapshot& snapshot) -> void
{
    rate_.store(snapshot.rate, std::memory_order_relaxed);
    eta_seconds_.store(snapshot.eta.count(), std::memory_order_relaxed);
    current_.store(snapshot.current, std::memory_order_relaxed);
}

auto Progress::set_prefix(std::string_view prefix) -> void
{
    const std::scoped_lock lock{prefix_mutex_};
    prefix_.assign(prefix);
}

auto Progress::finish() -> void
{
    if (!renderer_.joinable())
    {
        return;
    }
    renderer_.request_stop();
    renderer_.join();
}

auto Progress::format_line(std::string& output, std::string_view spinner) const
    -> void
{
    const std::size_t current = current_.load(std::memory_order_relaxed);
    const double rate = rate_.load(std::memory_order_relaxed);
    const std::int64_t eta_seconds
        = eta_seconds_.load(std::memory_order_relaxed);

    output.clear();
    output.append("\r\033[K");

    if (should_colorize())
    {
        fmt::format_to(
            std::back_inserter(output),
            style_for(ColorRole::success),
            "{}",
            spinner);
    }
    else
    {
        output.append(spinner);
    }
    output.push_back(' ');

    {
        const std::scoped_lock lock{prefix_mutex_};
        if (!prefix_.empty())
        {
            output.append(prefix_);
            output.push_back(' ');
        }
    }

    output.push_back('[');
    append_bar(output, current, total_);
    const double percentage
        = 100.0 * static_cast<double>(current) / static_cast<double>(total_);
    fmt::format_to(std::back_inserter(output), "] {:.1f}%", percentage);

    if (rate <= 0.0)
    {
        fmt::format_to(std::back_inserter(output), " | -- {}/s", rate_unit_);
    }
    else
    {
        fmt::format_to(
            std::back_inserter(output),
            " | {} {}/s",
            cli::AbbrNumber(rate),
            rate_unit_);
    }

    if (eta_seconds < 0)
    {
        output.append(" | ETA --:--:--");
    }
    else
    {
        const std::int64_t eta_hours = eta_seconds / 3600;
        const std::int64_t eta_minutes = eta_seconds % 3600 / 60;
        const std::int64_t eta_remainder_seconds = eta_seconds % 60;
        fmt::format_to(
            std::back_inserter(output),
            " | ETA {:02}:{:02}:{:02}",
            eta_hours,
            eta_minutes,
            eta_remainder_seconds);
    }
}

auto Progress::render(const std::stop_token& stop_token) -> void
{
    std::string frame;
    {
        const std::scoped_lock lock{prefix_mutex_};
        frame.reserve(prefix_.size() + 96);
    }

    while (!stop_token.stop_requested())
    {
        for (std::string_view spinner : spinner_frames)
        {
            if (stop_token.stop_requested())
            {
                break;
            }
            format_line(frame, spinner);
            std::fwrite(frame.data(), 1, frame.size(), stdout);
            std::fflush(stdout);
            std::this_thread::sleep_for(refresh_interval);
        }
    }

    format_line(frame, " ");
    frame.push_back('\n');
    std::fwrite(frame.data(), 1, frame.size(), stdout);
    std::fflush(stdout);
}

}  // namespace cli
