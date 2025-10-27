#include "gelex/utils/utils.h"

#include <chrono>

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <omp.h>
#include <Eigen/Dense>

namespace gelex
{

Timer::~Timer()
{
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start_;
    elapsed_time_ = duration.count();
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
