#include "logger_utils.h"

#include <iostream>
#include <string_view>

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
#include <spdlog/pattern_formatter.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace gelex
{
namespace detail
{

std::string with_std(double value, double std)
{
    return fmt::format("{:.6f} \u00B1 {:.4f}", value, std);
}

std::string title(const std::string& text, size_t total_length)
{
    if (text.empty())
    {
        return std::string(total_length, '=');  // Fallback to ASCII hyphen
    }

    size_t text_length = text.length();
    if (text_length >= total_length)
    {
        return text;
    }

    size_t remaining_space = total_length - text_length;
    size_t left_hyphens = remaining_space / 2;
    size_t right_hyphens = remaining_space - left_hyphens;  // Handles odd cases

    std::string left(left_hyphens, '=');
    std::string right(right_hyphens, '=');
    return left + text + right;
}

std::string subtitle(const std::string& str)
{
    return fmt::format(
        "{}", fmt::styled("[" + str + "]", fmt::fg(fmt::rgb(0x982756))));
}

std::string ToLowercase(std::string_view input)
{
    std::string result(input.begin(), input.end());
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}
std::string sigma_prior(std::string_view subscript, double nu, double s2)
{
    return fmt::format("σ²{} ∼ Inv-χ²(ν={:.2f}, S²={:.4f})", subscript, nu, s2);
}

std::string sigma_squared(std::string_view subscript)
{
    return fmt::format("σ²{}", subscript);
}

std::string h2(std::string_view subscript)
{
    return fmt::format("h²{}", subscript);
}

LevelFormatter::LevelFormatter()
{
    info_formatter_ = std::make_unique<spdlog::pattern_formatter>("%v");
    default_formatter_
        = std::make_unique<spdlog::pattern_formatter>("[%^%l%$] %v");
}

void LevelFormatter::format(
    const spdlog::details::log_msg& msg,
    spdlog::memory_buf_t& dest)
{
    if (msg.level == spdlog::level::info)
    {
        info_formatter_->format(msg, dest);
    }
    else
    {
        default_formatter_->format(msg, dest);
    }
}

std::unique_ptr<spdlog::formatter> LevelFormatter::clone() const
{
    return std::make_unique<LevelFormatter>();
}

Logger::Logger()
{
    try
    {
        auto console_sink
            = std::make_shared<spdlog::sinks::stdout_color_sink_st>();
        console_sink->set_level(spdlog::level::debug);
        console_sink->set_formatter(std::make_unique<LevelFormatter>());

        auto now = std::chrono::system_clock::now();
        auto time_t_now = std::chrono::system_clock::to_time_t(now);
        std::tm tm_now;
        localtime_r(&time_t_now, &tm_now);
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_st>(
            fmt::format("gelex_{:%F_%H-%M-%S}.log", tm_now));

        file_sink->set_level(spdlog::level::trace);
        file_sink->set_pattern("[%Y-%m-%d %H:%M:%S] [%l] %v");

        logger_ = std::make_shared<spdlog::logger>(
            "default", spdlog::sinks_init_list{console_sink, file_sink});
        spdlog::register_logger(logger_);
        logger_->set_level(spdlog::level::debug);
    }
    catch (const spdlog::spdlog_ex& ex)
    {
        std::cout << "Log initialization failed: " << ex.what() << "\n";
    }
}

std::shared_ptr<spdlog::logger> Logger::GetSpdLogger()
{
    return logger_;
}

std::shared_ptr<spdlog::logger> Logger::logger()
{
    static Logger instance;
    return instance.GetSpdLogger();
}

}  // namespace detail
}  // namespace gelex
