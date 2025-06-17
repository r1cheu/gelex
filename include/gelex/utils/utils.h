#pragma once

#include <chrono>
#include <memory>
#include <string>

#include <fmt/color.h>
#include <fmt/ranges.h>
#include <spdlog/common.h>
#include <spdlog/formatter.h>
#include <spdlog/logger.h>
#include <spdlog/spdlog.h>

namespace gelex
{

class LevelFormatter : public spdlog::formatter
{
   private:
    std::unique_ptr<spdlog::formatter> info_formatter_;
    std::unique_ptr<spdlog::formatter> default_formatter_;

   public:
    LevelFormatter();
    void format(const spdlog::details::log_msg& msg, spdlog::memory_buf_t& dest)
        override;
    std::unique_ptr<formatter> clone() const override;
};

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
