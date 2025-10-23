#pragma once

#include <atomic>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "Eigen/Core"
#include "gelex/barkeep.h"

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
        Eigen::Index n_iters,
        std::vector<std::atomic_ptrdiff_t>& progress_counters);

    template <typename T>
    void update(size_t chain_index, const std::string& status_name, T value)
    {
        if (chain_index >= num_chains_)
        {
            return;
        }

        // Directly update relevant fields for compact display
        if (status_name == "σ²_add" || status_name == "h²_add"
            || status_name == "σ²_dom" || status_name == "h²_dom"
            || status_name == "σ²_e" || status_name.rfind("π_", 0) == 0)
        {
            current_values_[status_name] = static_cast<double>(value);
            update_compact_status(chain_index);
        }
    }

    void show();
    void done();

   private:
    void update_compact_status(size_t chain_index);
    size_t num_chains_;

    std::vector<std::shared_ptr<bk::ProgressBarDisplay<std::atomic_ptrdiff_t>>>
        progress_bars_;
    std::vector<std::vector<std::shared_ptr<bk::StatusDisplay>>> statuses_;
    std::unordered_map<std::string, double> current_values_;
    std::shared_ptr<bk::CompositeDisplay> main_indicator_;
};

}  // namespace detail
}  // namespace gelex
