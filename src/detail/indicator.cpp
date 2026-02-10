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

#include "gelex/detail/indicator.h"
#include <atomic>
#include <vector>

#include <barkeep.h>
#include <fmt/format.h>
#include <Eigen/Core>

namespace gelex
{
namespace detail
{
namespace bk = barkeep;

auto create_progress_bar(
    std::atomic<size_t>& counter,
    size_t total,
    std::string_view format) -> ProgressBarDisplay
{
    std::vector<std::shared_ptr<bk::BaseDisplay>> elements;

    elements.push_back((bk::Animation(
        {.message = " ",
         .style = gelex::detail::GREEN_SPINNER,
         .interval = 0.08,
         .show = false})));
    auto before = bk::Status({.style = bk::Strings{" "}, .show = false});
    elements.push_back(before);
    elements.push_back(
        bk::ProgressBar(
            &counter,
            {.total = total,
             .format = std::string(format),
             .speed = 0.1,
             .style = gelex::detail::BAR_STYLE,
             .show = false}));
    auto after = bk::Status({.style = bk::Strings{" "}, .show = false});
    elements.push_back(after);

    return {
        .display = bk::Composite(elements, ""),
        .before = before,
        .after = after};
}

Indicator::Indicator(
    Eigen::Index n_iter,
    std::atomic_ptrdiff_t& progress_counter)
    : current_values_()
{
    std::vector<std::shared_ptr<bk::BaseDisplay>> displays;
    auto anim = bk::Animation(
        {.message = " ",
         .style = gelex::detail::GREEN_SPINNER,
         .interval = 0.08,
         .show = false});
    displays.emplace_back(anim);

    progress_bar_ = bk::ProgressBar(
        &progress_counter,
        {.total = n_iter,
         .format = "{bar} {value}/{total} [{speed:.1f}/s]",
         .speed = 0.1,
         .style = gelex::detail::BAR_STYLE,
         .show = false});
    displays.emplace_back(progress_bar_);

    status_ = bk::Status(
        {.message = "--", .style = bk::Strings{""}, .show = false});
    displays.emplace_back(status_);

    main_indicator_ = bk::Composite(displays, " ");
}

auto Indicator::update(StatusMetric metric, double value) -> void
{
    current_values_[status_metric_index(metric)] = value;
}

auto Indicator::flush_status() -> void
{
    update_compact_status();
}

auto Indicator::status_metric_index(StatusMetric metric) -> size_t
{
    return static_cast<size_t>(metric);
}

auto Indicator::format_status_line(const StatusSnapshot& values) -> std::string
{
    fmt::memory_buffer line_buffer;
    auto out_it = std::back_inserter(line_buffer);
    bool has_previous = false;

    auto append_metric = [&](std::string_view label, double value) -> void
    {
        if (has_previous)
        {
            out_it = fmt::format_to(out_it, " | ");
        }
        out_it = fmt::format_to(out_it, "{}: {:.3f}", label, value);
        has_previous = true;
    };

    if (const auto& add_h2
        = values[status_metric_index(StatusMetric::additive_heritability)])
    {
        append_metric("h²", *add_h2);
    }
    if (const auto& dom_h2
        = values[status_metric_index(StatusMetric::dominant_heritability)])
    {
        append_metric("δ²", *dom_h2);
    }
    if (const auto& res_var
        = values[status_metric_index(StatusMetric::residual_variance)])
    {
        append_metric("σ²_e", *res_var);
    }

    if (!has_previous)
    {
        return "--";
    }

    return fmt::to_string(line_buffer);
}

auto Indicator::update_compact_status() -> void
{
    status_->message(format_status_line(current_values_));
}

auto Indicator::show() -> void
{
    main_indicator_->show();
}

auto Indicator::done() -> void
{
    main_indicator_->done();
}

}  // namespace detail
}  // namespace gelex
