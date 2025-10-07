#include "indicator.h"
#include <fmt/format.h>
#include <vector>

#include "../../logger/logger_utils.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{
namespace detail
{
Indicator::Indicator(
    const BayesModel& model,
    size_t n_iter,
    std::vector<std::atomic<size_t>>& progress_counters)
    : num_chains_(progress_counters.size()),
      status_names_(create_status_names(model))
{
    std::vector<std::shared_ptr<bk::BaseDisplay>> all_chains_displays;

    for (int i = 0; i < status_names_.size(); ++i)
    {
        status_name_to_index_[status_names_[i]] = i;
    }

    statuses_.resize(num_chains_);

    for (int i = 0; i < num_chains_; ++i)
    {
        std::vector<std::shared_ptr<bk::BaseDisplay>> line_displays;

        auto anim = bk::Animation(
            {.style
             = bk::Strings{"⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"},
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

        statuses_[i].resize(status_names_.size());
        for (int j = 0; j < status_names_.size(); ++j)
        {
            auto status = bk::Status(
                {.message = status_names_[j] + ": --",
                 .style = bk::Strings{""},
                 .show = false});
            statuses_[i][j] = status;
            line_displays.push_back(status);
        }

        all_chains_displays.push_back(bk::Composite(line_displays, " "));
    }

    main_indicator_ = bk::Composite(all_chains_displays, "\n");
}

std::vector<std::string> Indicator::create_status_names(const BayesModel& model)
{
    std::vector<std::string> status_names;

    for (const auto& random : model.random())
    {
        status_names.emplace_back(sigma_squared("_" + random.name));
    }

    if (model.additive())
    {
        status_names.emplace_back(sigma_squared("_" + model.additive()->name));
        if (model.trait()->estimate_pi())
        {
            for (int i = 0; i < model.additive()->pi.size(); ++i)
            {
                status_names.emplace_back(fmt::format("π_{}", i));
            }
        }
    }

    status_names.emplace_back(sigma_squared("_e"));

    // Add heritability tracking for additive effects
    if (model.additive())
    {
        status_names.emplace_back(h2("_" + model.additive()->name));
    }

    return status_names;
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
