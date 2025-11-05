#pragma once

#include <atomic>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>

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
        std::span<std::atomic_ptrdiff_t> progress_counters);

    void
    update(size_t chain_index, const std::string& status_name, double value);

    void flush_status(size_t chain_index);

    void show();
    void done();

   private:
    void update_compact_status(size_t chain_index);
    size_t num_chains_;

    static const std::unordered_set<std::string_view> status_names_;

    std::vector<std::shared_ptr<bk::ProgressBarDisplay<std::atomic_ptrdiff_t>>>
        progress_bars_;
    std::vector<std::shared_ptr<bk::StatusDisplay>> statuses_;
    std::shared_ptr<bk::CompositeDisplay> main_indicator_;

    std::vector<std::map<std::string, double>> chain_values_;
    std::vector<std::atomic_bool> chain_dirty_flags_;
};

}  // namespace detail
}  // namespace gelex
