#pragma once

#include <chrono>
#include <memory>
#include <string>

#include <fmt/color.h>
#include <fmt/ranges.h>
#include <spdlog/logger.h>
#include <armadillo>

namespace gelex
{
using arma::dmat;
using arma::dvec;
using arma::sp_dmat;
template <typename Mat>
bool check_eye(const Mat& inputs)
{
    if (!inputs.is_square())
    {
        return false;
    }

    if (!inputs.is_diagmat())
    {
        return false;
    }

    for (size_t i = 0; i < inputs.n_rows; ++i)
    {
        if (inputs.at(i, i) != 1)
        {
            return false;
        }
    }
    return true;
}

std::string ToLowercase(std::string_view input);

// help functions for logging
std::string format_sigma_squared(const std::string& name);
std::string format_value_with_std(double value, double std);
std::string title(const std::string& text, size_t total_length = 80);
std::string subtitle(const std::string& str);
std::string item(const std::string& item);
std::string subitem(const std::string& item);

// fmt color might be removed
template <typename T>
auto green(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::rgb(0x00FF00)));
}

template <typename T>
auto pink(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::pink));
}

template <typename T>
auto cyan(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::cyan));
}

template <typename T>
auto rebecca_purple(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::rebecca_purple));
}
template <typename T>
auto rebecca_purple_vec(const T& value)
{
    return fmt::styled(
        fmt::join(value, ", "), fmt::fg(fmt::color::rebecca_purple));
}

template <typename T>
auto gray(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::gray));
}

template <typename T>
auto gold(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::gold));
}

template <typename T>
auto red(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::red));
}

template <typename T>
auto wine_red(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::rgb(0x982756)));
}

class Timer
{
   public:
    Timer(const Timer&) = default;
    Timer(Timer&&) = delete;
    Timer& operator=(const Timer&) = delete;
    Timer& operator=(Timer&&) = delete;
    explicit Timer(double& elapsed_time)
        : elapsed_time_{elapsed_time},
          start_{std::chrono::high_resolution_clock::now()} {};

    ~Timer();

   private:
    double& elapsed_time_;
    std::chrono::high_resolution_clock::time_point start_;
};

class Logger
{
   public:
    static std::shared_ptr<spdlog::logger> logger();

   private:
    Logger();
    std::shared_ptr<spdlog::logger> GetSpdLogger();
    std::shared_ptr<spdlog::logger> logger_;
};

std::string compute_time_left(
    const std::chrono::high_resolution_clock::time_point& start,
    size_t iter,
    size_t total_iter);

}  // namespace gelex
