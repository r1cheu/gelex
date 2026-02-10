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

#ifndef GELEX_DETAIL_INDICATOR_H_
#define GELEX_DETAIL_INDICATOR_H_

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string_view>

#include <Eigen/Core>

#include "barkeep.h"
#include "gelex/data/genotype_processor.h"

namespace gelex
{
namespace detail
{

inline const barkeep::BarParts BAR_STYLE{
    .left = "[",
    .right = "]",
    .fill = {"\033[1;36m━\033[0m"},
    .empty = {"-"}};

inline const barkeep::Strings GREEN_SPINNER{
    "\033[32m⠁\033[0m", "\033[32m⠁\033[0m", "\033[32m⠉\033[0m",
    "\033[32m⠙\033[0m", "\033[32m⠚\033[0m", "\033[32m⠒\033[0m",
    "\033[32m⠂\033[0m", "\033[32m⠂\033[0m", "\033[32m⠒\033[0m",
    "\033[32m⠲\033[0m", "\033[32m⠴\033[0m", "\033[32m⠤\033[0m",
    "\033[32m⠄\033[0m", "\033[32m⠄\033[0m", "\033[32m⠤\033[0m",
    "\033[32m⠠\033[0m", "\033[32m⠠\033[0m", "\033[32m⠤\033[0m",
    "\033[32m⠦\033[0m", "\033[32m⠖\033[0m", "\033[32m⠒\033[0m",
    "\033[32m⠐\033[0m", "\033[32m⠐\033[0m", "\033[32m⠒\033[0m",
    "\033[32m⠓\033[0m", "\033[32m⠋\033[0m", "\033[32m⠉\033[0m",
    "\033[32m⠈\033[0m", "\033[32m⠈\033[0m", "\033[32m \033[0m"};

struct ProgressBarDisplay
{
    std::shared_ptr<barkeep::CompositeDisplay> display;
    std::shared_ptr<barkeep::StatusDisplay> before;
    std::shared_ptr<barkeep::StatusDisplay> after;
};

auto create_progress_bar(
    std::atomic<size_t>& counter,
    size_t total,
    std::string_view format = "{bar}") -> ProgressBarDisplay;

class Indicator
{
   public:
    enum class StatusMetric : uint8_t
    {
        additive_heritability = 0,
        dominant_heritability = 1,
        residual_variance = 2,
    };

    Indicator(Eigen::Index n_iters, std::atomic_ptrdiff_t& progress_counter);

    auto update(StatusMetric metric, double value) -> void;

    auto flush_status() -> void;

    auto show() -> void;
    auto done() -> void;

   private:
    static constexpr size_t status_metric_count = 3;
    using StatusSnapshot
        = std::array<std::optional<double>, status_metric_count>;

    static auto status_metric_index(StatusMetric metric) -> size_t;
    static auto format_status_line(const StatusSnapshot& values) -> std::string;

    auto update_compact_status() -> void;

    std::shared_ptr<barkeep::ProgressBarDisplay<std::atomic_ptrdiff_t>>
        progress_bar_;
    std::shared_ptr<barkeep::StatusDisplay> status_;
    std::shared_ptr<barkeep::CompositeDisplay> main_indicator_;

    StatusSnapshot current_values_;
};

template <GenotypeProcessor Processor>
auto create_genotype_process_bar(int64_t& current, int64_t total)
    -> std::shared_ptr<barkeep::ProgressBarDisplay<int64_t>>
{
    return barkeep::ProgressBar(
        &current,
        {.total = total,
         .format = "      └─ {value}/{total} SNPs encoded "
                   "({speed:.1f} snp/s)",
         .speed = 0.1,
         .show = false});
};

}  // namespace detail
}  // namespace gelex

#endif  // GELEX_DETAIL_INDICATOR_H_
