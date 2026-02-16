/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
