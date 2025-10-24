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
    statuses_.reserve(num_chains_);

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

        auto compact_status = bk::Status(
            {.message = "--", .style = bk::Strings{""}, .show = false});
        statuses_.push_back(compact_status);
        line_displays.push_back(compact_status);

        all_chains_displays.push_back(bk::Composite(line_displays, " "));
    }

    main_indicator_ = bk::Composite(all_chains_displays, "\n");
}

void Indicator::update(
    size_t chain_index,
    const std::string& status_name,
    double value)
{
    if (chain_index >= num_chains_)
    {
        return;
    }

    if (status_name == "σ²_add" || status_name == "h²_add"
        || status_name == "σ²_dom" || status_name == "h²_dom"
        || status_name == "σ²_e" || status_name.rfind("π_", 0) == 0)
    {
        current_values_[status_name] = value;
        update_compact_status(chain_index);
    }
}

void Indicator::update_compact_status(size_t chain_index)
{
    if (chain_index >= num_chains_)
    {
        return;
    }

    std::string compact_line;

    auto add_var = current_values_.find("σ²_add");
    auto add_h2 = current_values_.find("h²_add");
    if (add_var != current_values_.end() && add_h2 != current_values_.end())
    {
        compact_line += fmt::format(
            "a(σ², h²): [{:.3f} {:.3f}] ", add_var->second, add_h2->second);
    }

    std::string pi_line = "π [";
    bool has_pi = false;
    for (size_t i = 0; i < 10; ++i)
    {
        auto it = current_values_.find(fmt::format("π_{}", i));
        if (it != current_values_.end())
        {
            if (has_pi)
                pi_line += ", ";
            pi_line += fmt::format("{:.3f}", it->second);
            has_pi = true;
        }
        else
        {
            break;
        }
    }
    if (has_pi)
    {
        compact_line += pi_line + "] ";
    }

    auto dom_var = current_values_.find("σ²_dom");
    auto dom_h2 = current_values_.find("h²_dom");
    if (dom_var != current_values_.end() && dom_h2 != current_values_.end())
    {
        compact_line += fmt::format(
            "d(σ², h²): [{:.3f} {:.3f}] ", dom_var->second, dom_h2->second);
    }

    auto res_var = current_values_.find("σ²_e");
    if (res_var != current_values_.end())
    {
        compact_line += fmt::format("σ²_e: {:.2f} ", res_var->second);
    }

    statuses_[chain_index]->message(compact_line);
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
