#ifndef GELEX_UTILS_UTILS_H_
#define GELEX_UTILS_UTILS_H_

#include <chrono>
#include <cstddef>  // 用于 size_t
#include <string>

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

class SmoothEtaCalculator
{
   public:
    explicit SmoothEtaCalculator(
        size_t total_items,
        double alpha = 0.1,
        size_t min_update_interval_ms = 500);

    std::string get_eta(size_t current_items);
    void reset(size_t new_total);
    std::string total_time_consumed() const;

   private:
    double calculate_eta_from_rate(size_t current, double rate) const;

    size_t total_items_;
    double alpha_;
    std::chrono::milliseconds min_update_interval_;

    std::chrono::steady_clock::time_point start_time_;
    std::chrono::steady_clock::time_point last_time_;
    size_t last_items_;

    double smooth_rate_;
    bool is_first_update_;
    double cached_eta_seconds_;
};

}  // namespace gelex

#endif  // GELEX_UTILS_UTILS_H_
