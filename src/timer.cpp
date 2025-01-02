#include "chenx/timer.h"

namespace chenx
{

Timer::~Timer()
{
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start_;
    elapsed_time_ = duration.count();
}
}  // namespace chenx
