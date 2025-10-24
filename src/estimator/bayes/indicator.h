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
    Indicator(
        Eigen::Index n_iters,
        std::vector<std::atomic_ptrdiff_t>& progress_counters);

    void
    update(size_t chain_index, const std::string& status_name, double value);

    void show();
    void done();

   private:
    void update_compact_status(size_t chain_index);
    size_t num_chains_;

    std::vector<std::shared_ptr<bk::ProgressBarDisplay<std::atomic_ptrdiff_t>>>
        progress_bars_;
    std::vector<std::shared_ptr<bk::StatusDisplay>> statuses_;
    std::unordered_map<std::string, double> current_values_;
    std::shared_ptr<bk::CompositeDisplay> main_indicator_;
};

}  // namespace detail
}  // namespace gelex
