#pragma once

#include <chrono>
#include <string>

#include <fmt/color.h>
#include <fmt/ranges.h>
#include <spdlog/common.h>
#include <spdlog/formatter.h>
#include <spdlog/logger.h>
#include <spdlog/spdlog.h>
#include <armadillo>

namespace gelex
{

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

std::string compute_time_left(
    const std::chrono::high_resolution_clock::time_point& start,
    size_t iter,
    size_t total_iter);

}  // namespace gelex
