#include "logger_utils.h"

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

}  // namespace detail
}  // namespace gelex
