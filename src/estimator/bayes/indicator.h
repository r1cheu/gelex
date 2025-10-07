#pragma once

#include <atomic>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "gelex/barkeep.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{
namespace detail
{
namespace bk = barkeep;
const bk::BarParts BAR_STYLE{
    .left = "[",
    .right = "]",
    .fill = {"\033[1;33m━\033[0m"},
    .empty = {"-"}};

class Indicator
{
   public:
    // Static progress bar style for reuse across the project
    Indicator(
        const BayesModel& model,
        size_t n_iters,
        std::vector<std::atomic<size_t>>& progress_counters);

    template <typename T>
    void update(size_t chain_index, const std::string& status_name, T value)
    {
        if (chain_index >= num_chains_)
        {
            return;
        }
        auto it = status_name_to_index_.find(status_name);
        if (it == status_name_to_index_.end())
        {
            return;
        }
        size_t status_idx = it->second;

        statuses_[chain_index][status_idx]->message(
            fmt::format("{}: {:.4f}", status_name, value));
    }

    void show();
    void done();

   private:
    static std::vector<std::string> create_status_names(
        const BayesModel& model);

    size_t num_chains_;
    std::vector<std::string> status_names_;
    std::map<std::string, size_t> status_name_to_index_;

    std::vector<std::shared_ptr<bk::ProgressBarDisplay<std::atomic_size_t>>>
        progress_bars_;
    std::vector<std::vector<std::shared_ptr<bk::StatusDisplay>>> statuses_;
    std::shared_ptr<bk::CompositeDisplay> main_indicator_;
};

}  // namespace detail
}  // namespace gelex
