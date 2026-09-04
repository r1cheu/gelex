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

#ifndef APPS_CLI_PROGRESS_H_
#define APPS_CLI_PROGRESS_H_

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
#include <string_view>
#include <thread>

namespace cli
{

inline auto make_rate(std::size_t initial_current = 0)
{
    return [last_time = std::chrono::steady_clock::now(),
            last_current = initial_current,
            smoothed_rate = 0.0,
            has_rate = false](std::size_t current) mutable
    {
        constexpr double smoothing_factor = 0.2;
        const auto now = std::chrono::steady_clock::now();
        const std::size_t completed = current - last_current;
        const double elapsed_seconds
            = std::chrono::duration<double>(now - last_time).count();
        const double instantaneous_rate
            = static_cast<double>(completed) / elapsed_seconds;
        smoothed_rate = has_rate
                            ? (smoothing_factor * instantaneous_rate)
                                  + ((1.0 - smoothing_factor) * smoothed_rate)
                            : instantaneous_rate;
        last_time = now;
        last_current = current;
        has_rate = true;
        return smoothed_rate;
    };
}

inline auto make_eta(std::size_t total)
{
    return [total,
            start_time = std::chrono::steady_clock::now()](std::size_t current)
    {
        if (current == total)
        {
            return std::chrono::seconds{0};
        }

        const auto now = std::chrono::steady_clock::now();
        const double elapsed_seconds
            = std::chrono::duration<double>(now - start_time).count();
        const double average_rate
            = static_cast<double>(current) / elapsed_seconds;
        const auto eta = std::chrono::duration<double>{
            static_cast<double>(total - current) / average_rate};
        return std::chrono::ceil<std::chrono::seconds>(eta);
    };
}

struct ProgressSnapshot
{
    std::size_t current;
    double rate;
    std::chrono::seconds eta;
};

/**
 * @brief Renders progress on a single terminal line.
 *
 * A background thread renders the latest snapshot when standard output is a
 * TTY and `GELEX_NO_TTY` is not set. The owning thread publishes progress with
 * update() and stops rendering with finish() or destruction.
 */
class Progress
{
   public:
    /**
     * @brief Constructs a progress display and starts rendering when enabled.
     *
     * @param prefix Text displayed before the progress bar.
     * @param total Final progress value.
     * @param rate_unit Unit displayed for the rate per second.
     */
    Progress(std::string prefix, std::size_t total, std::string rate_unit);

    Progress(const Progress&) = delete;
    Progress(Progress&&) = delete;
    auto operator=(const Progress&) -> Progress& = delete;
    auto operator=(Progress&&) -> Progress& = delete;

    ~Progress();

    /** @brief Publishes the latest progress values to the renderer. */
    auto update(const ProgressSnapshot& snapshot) -> void;

    /** @brief Replaces the text displayed before the progress bar. */
    auto set_prefix(std::string_view prefix) -> void;

    /** @brief Stops rendering and writes the final line when enabled. */
    auto finish() -> void;

   private:
    auto format_line(std::string& output, std::string_view spinner) const
        -> void;
    auto render(const std::stop_token& stop_token) -> void;

    std::string prefix_;
    mutable std::mutex prefix_mutex_;
    std::string rate_unit_;
    std::size_t total_;
    std::atomic_size_t current_{0};
    std::atomic<double> rate_{0.0};
    std::atomic<std::int64_t> eta_seconds_{-1};
    std::jthread renderer_;
};

}  // namespace cli

#endif  // APPS_CLI_PROGRESS_H_
