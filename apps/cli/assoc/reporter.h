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

#ifndef APPS_CLI_ASSOC_REPORTER_H_
#define APPS_CLI_ASSOC_REPORTER_H_

#include <Eigen/Core>
#include <cstddef>
#include <string_view>
#include <vector>

#include "gelex/freq/reml/summary.h"

#include "cli/progress.h"

namespace gelex
{
class FreqModel;
}  // namespace gelex

namespace cli
{

class AssocReporter
{
   public:
    explicit AssocReporter(std::size_t total_snps);

    static auto show_dataset_summary(
        const gelex::FreqModel& model,
        Eigen::Index n_snps) -> void;
    static auto show_scan_header(size_t total_snps, int chunk_size, bool loco)
        -> void;
    auto update_scan_progress(size_t current, double rate) -> void;
    auto show_loco_phase(std::string_view phase) -> void;
    auto show_loco_reml_summary(
        const std::vector<gelex::LocoRemlResult>& results) -> void;
    auto finish_scan() -> void;

   private:
    cli::Progress progress_;
    decltype(cli::make_eta(std::size_t{})) estimate_eta_;
};

}  // namespace cli

#endif  // APPS_CLI_ASSOC_REPORTER_H_
