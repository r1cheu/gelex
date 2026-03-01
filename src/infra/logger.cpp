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

#include "gelex/infra/logger.h"

#include <memory>

#include <spdlog/pattern_formatter.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace
{
std::shared_ptr<spdlog::logger> g_logger = nullptr;
std::once_flag g_init_flag;

class LevelFormatter : public spdlog::formatter
{
   public:
    LevelFormatter()
    {
        static constexpr std::string_view info_pattern = "%v";
        static constexpr std::string_view default_pattern = "[%^%l%$] %v";

        info_formatter_ = std::make_unique<spdlog::pattern_formatter>(
            std::string(info_pattern));
        default_formatter_ = std::make_unique<spdlog::pattern_formatter>(
            std::string(default_pattern));
    }

    void format(const spdlog::details::log_msg& msg, spdlog::memory_buf_t& dest)
        override
    {
        if (msg.level != spdlog::level::info)
        {
            default_formatter_->format(msg, dest);
        }
        else
        {
            info_formatter_->format(msg, dest);
        }
    }

    [[nodiscard]] std::unique_ptr<spdlog::formatter> clone() const override
    {
        return std::make_unique<LevelFormatter>();
    }

   private:
    std::unique_ptr<spdlog::formatter> info_formatter_;
    std::unique_ptr<spdlog::formatter> default_formatter_;
};

// Strips ANSI CSI escape sequences (e.g. fmt::fg color codes) and bare CR
// from msg.payload before formatting, so log files remain clean and readable.
class StripAnsiFormatter : public spdlog::formatter
{
   public:
    explicit StripAnsiFormatter(std::string pattern)
        : pattern_(std::move(pattern)),
          formatter_(std::make_unique<spdlog::pattern_formatter>(pattern_))
    {
    }

    void format(const spdlog::details::log_msg& msg, spdlog::memory_buf_t& dest)
        override
    {
        std::string clean;
        clean.reserve(msg.payload.size());
        for (std::size_t i = 0; i < msg.payload.size();)
        {
            // ESC + '[' marks the start of a CSI sequence; skip until final
            // byte (0x40–0x7E per ECMA-48)
            if (msg.payload[i] == '\x1b' && i + 1 < msg.payload.size()
                && msg.payload[i + 1] == '[')
            {
                i += 2;
                while (i < msg.payload.size()
                       && (msg.payload[i] < 0x40 || msg.payload[i] > 0x7e))
                {
                    ++i;
                }
                if (i < msg.payload.size())
                {
                    ++i;
                }
            }
            else if (msg.payload[i] == '\r')
            {
                ++i;
            }
            else
            {
                clean.push_back(msg.payload[i++]);
            }
        }

        spdlog::details::log_msg clean_msg(msg);
        clean_msg.payload = spdlog::string_view_t(clean.data(), clean.size());
        formatter_->format(clean_msg, dest);
    }

    [[nodiscard]] std::unique_ptr<spdlog::formatter> clone() const override
    {
        return std::make_unique<StripAnsiFormatter>(pattern_);
    }

   private:
    std::string pattern_;
    std::unique_ptr<spdlog::pattern_formatter> formatter_;
};
}  // namespace

namespace gelex::logging
{

void initialize(std::string_view output_prefix)
{
    std::call_once(
        g_init_flag,
        [&]()
        {
            auto console_sink
                = std::make_shared<spdlog::sinks::stdout_color_sink_st>();
            console_sink->set_level(spdlog::level::debug);
            console_sink->set_formatter(std::make_unique<LevelFormatter>());

            const std::string log_filename
                = std::format("{}.log", output_prefix);
            auto file_sink
                = std::make_shared<spdlog::sinks::basic_file_sink_st>(
                    log_filename, true);
            file_sink->set_level(spdlog::level::trace);
            file_sink->set_formatter(
                std::make_unique<StripAnsiFormatter>(
                    "[%Y-%m-%d %H:%M:%S.%e] [%l] %v"));

            g_logger = std::make_shared<spdlog::logger>(
                "gelex_logger",
                spdlog::sinks_init_list{console_sink, file_sink});

            g_logger->set_level(spdlog::level::trace);
            g_logger->flush_on(spdlog::level::err);

            spdlog::register_logger(g_logger);
        });
}

std::shared_ptr<spdlog::logger>& get()
{
    return g_logger;
}

}  // namespace gelex::logging
