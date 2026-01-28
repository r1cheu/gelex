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

#ifndef GELEX_ESTIMATOR_BAYES_INDICATOR_H_
#define GELEX_ESTIMATOR_BAYES_INDICATOR_H_

#include <atomic>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>

#include "barkeep.h"
#include "gelex/data/variant_processor.h"

namespace gelex
{
namespace detail
{
const barkeep::BarParts BAR_STYLE{
    .left = "[",
    .right = "]",
    .fill = {"\033[1;36m━\033[0m"},
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

    std::vector<
        std::shared_ptr<barkeep::ProgressBarDisplay<std::atomic_ptrdiff_t>>>
        progress_bars_;
    std::vector<std::shared_ptr<barkeep::StatusDisplay>> statuses_;
    std::shared_ptr<barkeep::CompositeDisplay> main_indicator_;

    std::vector<std::map<std::string, double>> chain_values_;
    std::vector<std::atomic_bool> chain_dirty_flags_;
};

template <VariantProcessor Processor>
auto create_genotype_process_bar(int64_t& current, int64_t total)
    -> std::shared_ptr<barkeep::ProgressBarDisplay<int64_t>>
{
    return barkeep::ProgressBar(
        &current,
        {.total = total,
         .format = "      └─ {value}/{total} SNPs encoded "
                   "({speed:.1f} snp/s)",
         .speed = 0.1,
         .show = false});
};

}  // namespace detail
}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_INDICATOR_H_
