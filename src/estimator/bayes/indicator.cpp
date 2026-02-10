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

#include "indicator.h"
#include <atomic>
#include <optional>
#include <vector>

#include <barkeep.h>
#include <fmt/format.h>
#include <Eigen/Core>

namespace gelex
{
namespace detail
{
namespace bk = barkeep;

const std::unordered_set<std::string_view> Indicator::status_names_{
    "additive_heritability",
    "dominant_heritability",
    "residual_variance"};

Indicator::Indicator(
    Eigen::Index n_iter,
    std::atomic_ptrdiff_t& progress_counter)
    : dirty_flag_(false)
{
    std::vector<std::shared_ptr<bk::BaseDisplay>> displays;
    auto anim = bk::Animation(
        {.message = " ",
         .style = bk::Strings{"\033[32m⠁\033[0m", "\033[32m⠁\033[0m",
                              "\033[32m⠉\033[0m", "\033[32m⠙\033[0m",
                              "\033[32m⠚\033[0m", "\033[32m⠒\033[0m",
                              "\033[32m⠂\033[0m", "\033[32m⠂\033[0m",
                              "\033[32m⠒\033[0m", "\033[32m⠲\033[0m",
                              "\033[32m⠴\033[0m", "\033[32m⠤\033[0m",
                              "\033[32m⠄\033[0m", "\033[32m⠄\033[0m",
                              "\033[32m⠤\033[0m", "\033[32m⠠\033[0m",
                              "\033[32m⠠\033[0m", "\033[32m⠤\033[0m",
                              "\033[32m⠦\033[0m", "\033[32m⠖\033[0m",
                              "\033[32m⠒\033[0m", "\033[32m⠐\033[0m",
                              "\033[32m⠐\033[0m", "\033[32m⠒\033[0m",
                              "\033[32m⠓\033[0m", "\033[32m⠋\033[0m",
                              "\033[32m⠉\033[0m", "\033[32m⠈\033[0m",
                              "\033[32m⠈\033[0m", "\033[32m \033[0m"},
         .interval = 0.08,
         .show = false});
    displays.emplace_back(anim);

    progress_bar_ = bk::ProgressBar(
        &progress_counter,
        {.total = n_iter,
         .format = "{bar} {value}/{total} [{speed:.1f}/s]",
         .speed = 0.1,
         .style = BAR_STYLE,
         .show = false});
    displays.emplace_back(progress_bar_);

    status_ = bk::Status(
        {.message = "--", .style = bk::Strings{""}, .show = false});
    displays.emplace_back(status_);

    main_indicator_ = bk::Composite(displays, " ");
}

void Indicator::update(const std::string& status_name, double value)
{
    std::string_view name_view = status_name;
    if (status_names_.contains(name_view))
    {
        current_values_[status_name] = value;
        dirty_flag_.store(true, std::memory_order_relaxed);
    }
}

void Indicator::flush_status()
{
    if (dirty_flag_.exchange(false, std::memory_order_acquire))
    {
        update_compact_status();
    }
}

void Indicator::update_compact_status()
{
    fmt::memory_buffer line_buffer;
    auto out_it = std::back_inserter(line_buffer);

    auto get_val = [&](const std::string& key) -> std::optional<double>
    {
        auto it = current_values_.find(key);
        if (it != current_values_.end())
        {
            return it->second;
        }
        return std::nullopt;
    };

    if (auto add_h2 = get_val("additive_heritability"))
    {
        out_it = fmt::format_to(out_it, "h²: {:.3f} ", *add_h2);
    }
    if (auto dom_h2 = get_val("dominant_heritability"))
    {
        out_it = fmt::format_to(out_it, " δ²: {:.3f}", *dom_h2);
    }
    if (auto res_var = get_val("residual_variance"))
    {
        out_it = fmt::format_to(out_it, " σ²_e: {:.3f} ", *res_var);
    }

    status_->message(fmt::to_string(line_buffer));
}

void Indicator::show()
{
    main_indicator_->show();
}
void Indicator::done()
{
    main_indicator_->done();
}

}  // namespace detail
}  // namespace gelex
