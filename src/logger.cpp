#include "gelex/logger.h"

#include <memory>

#include <spdlog/pattern_formatter.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace
{
// Use an anonymous namespace to hide global state
std::shared_ptr<spdlog::logger> g_logger = nullptr;
std::once_flag g_init_flag;
}  // namespace

namespace gelex::logging
{

LevelFormatter::LevelFormatter()
{
    static constexpr std::string_view info_pattern = "%v";
    static constexpr std::string_view default_pattern = "[%^%l%$] %v";

    info_formatter_ = std::make_unique<spdlog::pattern_formatter>(
        std::string(info_pattern));
    default_formatter_ = std::make_unique<spdlog::pattern_formatter>(
        std::string(default_pattern));
}

void LevelFormatter::format(
    const spdlog::details::log_msg& msg,
    spdlog::memory_buf_t& dest)
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

std::unique_ptr<spdlog::formatter> LevelFormatter::clone() const
{
    return std::make_unique<LevelFormatter>();
}

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
            file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");

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
