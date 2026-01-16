#include "utils.h"

#include <chrono>

namespace gelex
{

Timer::~Timer()
{
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start_;
    elapsed_time_ = duration.count();
}
}  // namespace gelex
