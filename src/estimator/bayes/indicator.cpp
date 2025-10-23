#include "indicator.h"
#include <fmt/format.h>
#include <atomic>
#include <vector>
#include "Eigen/Core"

namespace gelex
{
namespace detail
{
Indicator::Indicator(
    Eigen::Index n_iter,
    std::vector<std::atomic_ptrdiff_t>& progress_counters)
    : num_chains_(progress_counters.size())
{
    std::vector<std::shared_ptr<bk::BaseDisplay>> all_chains_displays;

    statuses_.resize(num_chains_);

    for (size_t i = 0; i < num_chains_; ++i)
    {
        std::vector<std::shared_ptr<bk::BaseDisplay>> line_displays;

        auto anim = bk::Animation(
            {.style = bk::Strings{"⣷", "⣯", "⣟", "⡿", "⢿", "⣻", "⣽", "⣾"},
             .interval = 0.08,
             .show = false});
        line_displays.push_back(anim);

        auto pbar = bk::ProgressBar(
            &progress_counters[i],
            {.total = n_iter,
             .format = fmt::format("{}", i + 1)
                       + " {bar} {value}/{total} ({speed:.1f}/s)",
             .speed = 0.1,
             .style = BAR_STYLE,
             .show = false});
        progress_bars_.push_back(pbar);
        line_displays.push_back(pbar);

        // Use a single compact status display instead of multiple ones
        auto compact_status = bk::Status(
            {.message = "--", .style = bk::Strings{""}, .show = false});
        statuses_[i].push_back(compact_status);
        line_displays.push_back(compact_status);

        all_chains_displays.push_back(bk::Composite(line_displays, " "));
    }

    main_indicator_ = bk::Composite(all_chains_displays, "\n");
}

void Indicator::update_compact_status(size_t chain_index)
{
    if (chain_index >= num_chains_ || statuses_[chain_index].empty())
    {
        return;
    }

    // Build compact single line display
    std::string compact_line;

    // Add additive variance and heritability
    if (current_values_.contains("σ²_additive")
        && current_values_.contains("h²_additive"))
    {
        compact_line += fmt::format(
            "a(σ², h²): [{:.3f} {:.3f}] ",
            current_values_.at("σ²_additive"),
            current_values_.at("h²_additive"));
    }

    // Add mixture parameters
    std::vector<double> pi_values;

    // Collect π values in order (π_0, π_1, π_2, ...)
    size_t i = 0;
    while (true)
    {
        std::string pi_key = fmt::format("π_{}", i);
        if (current_values_.contains(pi_key))
        {
            pi_values.push_back(current_values_.at(pi_key));
            ++i;
        }
        else
        {
            break;
        }
    }

    compact_line += "π [";
    for (size_t i = 0; i < pi_values.size(); ++i)
    {
        if (i > 0)
        {
            compact_line += ", ";
        }
        compact_line += fmt::format("{:.3f}", pi_values[i]);
    }
    compact_line += "] ";

    // Add dominance variance and heritability
    if (current_values_.contains("σ²_dominant")
        && current_values_.contains("h²_dominant"))
    {
        compact_line += fmt::format(
            "d(σ², h²): [{:.3f} {:.3f}] ",
            current_values_.at("σ²_dominant"),
            current_values_.at("h²_dominant"));
    }

    // Add residual variance
    if (current_values_.contains("σ²_e"))
    {
        compact_line
            += fmt::format("σ²_e: {:.2f} ", current_values_.at("σ²_e"));
    }

    // Update the compact status display
    statuses_[chain_index][0]->message(compact_line);
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
