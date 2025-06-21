#include "gelex/utils/formatter.h"

#include <string>
#include <string_view>

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <armadillo>

namespace gelex
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
    return fmt::format("{}", wine_red("[" + str + "]"));
}

std::string ToLowercase(std::string_view input)
{
    std::string result(input.begin(), input.end());
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}
std::string scale_inv_chisq(double nu, double s2)
{
    return fmt::format("Inv-\u03C7\u00B2(\u03BD={}, s\u00B2={:.4f})", nu, s2);
}

std::string sigma_squared(const std::string& subscript)
{
    return fmt::format("\u03C3\u00B2{}", subscript);
}

std::string h2(const std::string& subscript)
{
    return fmt::format("h\u00B2{}", subscript);
}

std::string sigma_prior(const std::string& subscript, double nu, double s2)
{
    return fmt::format(
        "{} \u223C {}", sigma_squared(subscript), scale_inv_chisq(nu, s2));
}
}  // namespace gelex
