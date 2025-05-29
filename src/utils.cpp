#include "gelex/utils.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <spdlog/logger.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <armadillo>

namespace gelex
{
std::string format_sigma_squared(const std::string& name)
{
    return fmt::format("\u03C3\u00B2_{}", name);
}

std::string format_value_with_std(double value, double std)
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
    return fmt::format("{}", wine_red("[" + str + "]"));
}

std::string item(const std::string& item)
{
    return fmt::format(" \u25AA {}", item);
}

std::string subitem(const std::string& item)
{
    return fmt::format(" - {}", item);
}

std::string ToLowercase(std::string_view input)
{
    std::string result(input.begin(), input.end());
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

Timer::~Timer()
{
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start_;
    elapsed_time_ = duration.count();
}

Logger::Logger()
{
    try
    {
        auto console_sink
            = std::make_shared<spdlog::sinks::stdout_color_sink_st>();
        console_sink->set_level(spdlog::level::debug);
        console_sink->set_pattern("[%^%l%$] %v");

        auto now = std::chrono::system_clock::now();
        auto file_sink
            = std::make_shared<spdlog::sinks::basic_file_sink_st>(fmt::format(
                "gelex_{:%F_%H-%M-%S}.log",
                fmt::localtime(std::chrono::system_clock::to_time_t(now))));

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

std::string compute_time_left(
    const std::chrono::high_resolution_clock::time_point& start,
    size_t iter,
    size_t total_iter)
{
    // Handle edge cases
    if (iter == 0)
    {
        return "--:--:--";
    }

    if (iter >= total_iter)
    {
        return "00:00:00";
    }

    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    auto time_per_iter = elapsed / iter;
    auto time_left = time_per_iter * (total_iter - iter);

    auto seconds_left
        = std::chrono::duration_cast<std::chrono::seconds>(time_left);
    return fmt::format("{:%H:%M:%S}", seconds_left);
}
}  // namespace gelex
