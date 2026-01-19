#include "indicator.h"
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

const std::unordered_set<std::string_view> Indicator::status_names_{
    "additive_heritability",
    "dominant_heritability",
    "residual_variance"};

Indicator::Indicator(
    Eigen::Index n_iter,
    std::span<std::atomic_ptrdiff_t> progress_counters)
    : num_chains_(progress_counters.size()),
      chain_values_(progress_counters.size()),
      chain_dirty_flags_(progress_counters.size())
{
    std::vector<std::shared_ptr<bk::BaseDisplay>> all_chains_displays;
    all_chains_displays.reserve(num_chains_);
    statuses_.reserve(num_chains_);
    progress_bars_.reserve(num_chains_);

    for (size_t i = 0; i < num_chains_; ++i)
    {
        chain_dirty_flags_[i].store(false, std::memory_order_relaxed);
    }

    for (size_t i = 0; i < num_chains_; ++i)
    {
        std::vector<std::shared_ptr<bk::BaseDisplay>> line_displays;

        auto anim = bk::Animation(
            {.message = " ",
             .style
             = bk::Strings{"⠁", "⠁", "⠉", "⠙", "⠚", "⠒", "⠂", "⠂", "⠒", "⠲",
                           "⠴", "⠤", "⠄", "⠄", "⠤", "⠠", "⠠", "⠤", "⠦", "⠖",
                           "⠒", "⠐", "⠐", "⠒", "⠓", "⠋", "⠉", "⠈", "⠈", " "},
             .interval = 0.08,
             .show = false});
        line_displays.emplace_back(anim);
        auto pbar = bk::ProgressBar(
            &progress_counters[i],
            {.total = n_iter,
             .format = "{bar} {value}/{total} [{speed:.1f}/s]",
             .speed = 0.1,
             .style = BAR_STYLE,
             .show = false});
        progress_bars_.emplace_back(pbar);
        line_displays.emplace_back(pbar);

        auto compact_status = bk::Status(
            {.message = "--", .style = bk::Strings{""}, .show = false});
        statuses_.emplace_back(compact_status);
        line_displays.emplace_back(compact_status);

        all_chains_displays.emplace_back(bk::Composite(line_displays, " "));
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

    std::string_view name_view = status_name;
    if (status_names_.contains(name_view))
    {
        chain_values_[chain_index][status_name] = value;
        chain_dirty_flags_[chain_index].store(true, std::memory_order_relaxed);
    }
}

void Indicator::flush_status(size_t chain_index)
{
    if (chain_index >= num_chains_)
    {
        return;
    }

    if (chain_dirty_flags_[chain_index].exchange(
            false, std::memory_order_acquire))
    {
        update_compact_status(chain_index);
    }
}

void Indicator::update_compact_status(size_t chain_index)
{
    fmt::memory_buffer line_buffer;
    auto out_it = std::back_inserter(line_buffer);

    const auto& current_values = chain_values_[chain_index];

    auto get_val = [&](const std::string& key) -> std::optional<double>
    {
        auto it = current_values.find(key);
        if (it != current_values.end())
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

    statuses_[chain_index]->message(fmt::to_string(line_buffer));
    if (chain_index >= num_chains_)
    {
        return;
    }
}

void Indicator::show()
{
    main_indicator_->show();
}
void Indicator::done()
{
    main_indicator_->done();
}

auto create_association_progress_bar(size_t& current, size_t total)
    -> AssociationPbar
{
    std::vector<std::shared_ptr<bk::BaseDisplay>> elements{bk::ProgressBar(
        &current,
        {.total = total,
         .format = "   {bar}",
         .style = BAR_STYLE,
         .show = false})};
    auto status = bk::Status(
        {.message = "--%  --:--:--", .style = bk::Strings{""}, .show = false});
    elements.push_back(status);
    return {.pbar = bk::Composite(elements, " "), .status = status};
};

}  // namespace detail
}  // namespace gelex
