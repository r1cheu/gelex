#pragma once

#include <chrono>

namespace gelex
{

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

}  // namespace gelex
