#include "utils.h"

#include <chrono>
#include "utils/formatter.h"

namespace gelex
{

Timer::~Timer()
{
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start_;
    elapsed_time_ = duration.count();
}

SmoothEtaCalculator::SmoothEtaCalculator(
    size_t total_items,
    double alpha,
    size_t min_update_interval_ms)
    : total_items_(total_items),
      alpha_(alpha),
      min_update_interval_(min_update_interval_ms),
      start_time_(std::chrono::steady_clock::now()),
      last_time_(start_time_),
      last_items_(0),
      smooth_rate_(0.0),
      is_first_update_(true),
      cached_eta_seconds_(0.0)
{
}

std::string SmoothEtaCalculator::get_eta(size_t current_items)
{
    auto now = std::chrono::steady_clock::now();

    auto elapsed_since_last
        = std::chrono::duration_cast<std::chrono::milliseconds>(
              now - last_time_)
              .count();

    if (elapsed_since_last < min_update_interval_.count() && !is_first_update_)
    {
        return gelex::format_eta(
            calculate_eta_from_rate(current_items, smooth_rate_));
    }

    std::chrono::duration<double> delta_time = now - last_time_;
    double delta_sec = delta_time.count();

    if (current_items < last_items_)
    {
        last_items_ = current_items;
        return gelex::format_eta(
            calculate_eta_from_rate(current_items, smooth_rate_));
    }

    size_t delta_items = current_items - last_items_;

    double instant_rate = 0.0;
    if (delta_sec > 1e-9)
    {
        instant_rate = static_cast<double>(delta_items) / delta_sec;
    }

    if (is_first_update_)
    {
        if (delta_sec > 0)
        {
            smooth_rate_ = instant_rate;
        }
        is_first_update_ = false;
    }
    else
    {
        smooth_rate_
            = (instant_rate * alpha_) + (smooth_rate_ * (1.0 - alpha_));
    }

    last_time_ = now;
    last_items_ = current_items;

    cached_eta_seconds_ = calculate_eta_from_rate(current_items, smooth_rate_);
    return gelex::format_eta(cached_eta_seconds_);
}
std::string SmoothEtaCalculator::total_time_consumed() const
{
    auto now = std::chrono::steady_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::duration<double>>(
                          now - start_time_)
                          .count();
    return gelex::format_eta(duration);
}

void SmoothEtaCalculator::reset(size_t new_total)
{
    total_items_ = new_total;
    start_time_ = std::chrono::steady_clock::now();
    last_time_ = start_time_;
    last_items_ = 0;
    smooth_rate_ = 0.0;
    is_first_update_ = true;
    cached_eta_seconds_ = 0.0;
}

double SmoothEtaCalculator::calculate_eta_from_rate(size_t current, double rate)
    const
{
    if (rate <= 1e-6)
    {
        return 999999.0;
    }
    if (current >= total_items_)
    {
        return 0.0;
    }
    return static_cast<double>(total_items_ - current) / rate;
}
}  // namespace gelex
