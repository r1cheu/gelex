#include "gelex/estimator/bayes/indicator.h"
#include <fmt/format.h>
#include <vector>
#include "gelex/model/bayes/policy.h"
#include "gelex/utils/formatter.h"
namespace gelex
{
Indicator::Indicator(
    size_t n_chains,
    size_t n_iters,
    std::vector<std::atomic<size_t>>& progress_counters,
    const std::vector<std::string>& status_names)
    : num_chains_(n_chains), status_names_(status_names)
{
    std::vector<std::shared_ptr<bk::BaseDisplay>> all_chains_displays;

    for (size_t i = 0; i < status_names.size(); ++i)
    {
        status_name_to_index_[status_names[i]] = i;
    }

    statuses_.resize(n_chains);

    for (size_t i = 0; i < n_chains; ++i)
    {
        std::vector<std::shared_ptr<bk::BaseDisplay>> line_displays;

        auto anim = bk::Animation(
            {.style
             = bk::Strings{"⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"},
             .interval = 0.08,
             .show = false});
        line_displays.push_back(anim);

        bk::BarParts custom_bar_style;
        custom_bar_style.left = "[";
        custom_bar_style.right = "]";
        custom_bar_style.fill = {"\033[1;33m━\033[0m"};
        custom_bar_style.empty = {"─"};

        auto pbar = bk::ProgressBar(
            &progress_counters[i],
            {.total = n_iters,
             .format = fmt::format("{}", i + 1)
                       + " {bar} {value}/{total} ({speed:.1f}/s)",
             .speed = 0.1,
             .style = custom_bar_style,
             .show = false});
        progress_bars_.push_back(pbar);
        line_displays.push_back(pbar);

        statuses_[i].resize(status_names.size());
        for (size_t j = 0; j < status_names.size(); ++j)
        {
            auto status = bk::Status(
                {.message = status_names[j] + ": --",
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

    for (const auto& genetic : model.genetic())
    {
        status_names.emplace_back(sigma_squared("_" + genetic.name));
        if (bayes_trait_estimate_pi[to_index(genetic.type)])
        {
            for (size_t i = 0; i < genetic.pi.size(); ++i)
            {
                status_names.emplace_back(
                    fmt::format("π{}_{}", i, genetic.name));
            }
        }
    }

    status_names.emplace_back(sigma_squared("_e"));

    for (const auto& genetic : model.genetic())
    {
        status_names.emplace_back(h2("_" + genetic.name));
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

}  // namespace gelex
