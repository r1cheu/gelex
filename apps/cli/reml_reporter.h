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

#ifndef APPS_CLI_REML_REPORTER_H_
#define APPS_CLI_REML_REPORTER_H_

#include <cstddef>
#include <string_view>
#include <variant>
#include <vector>

#include "gelex/algo/reml/summary.h"
#include "gelex/infra/logging/reml_event.h"

#include "cli/table.h"

namespace gelex
{
class FreqModel;
}  // namespace gelex

namespace cli
{

class RemlReporter
{
   public:
    auto show_dataset_summary(
        const gelex::FreqModel& model,
        std::string_view pheno_name) const -> void;
    auto on_event(const gelex::RemlIterationEvent& e) -> void;
    auto on_event(const gelex::RemlConstrainedEvent& e) -> void;
    auto show_result(
        const gelex::FreqModel& model,
        const gelex::RemlSummary& summary,
        size_t max_iter) const -> void;

    auto as_observer() -> gelex::RemlObserver
    {
        return [this](const gelex::RemlEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    Table iter_table_;
    bool header_printed_ = false;
};

void print_loco_reml_summary(const std::vector<gelex::LocoRemlResult>& results);

}  // namespace cli

#endif  // APPS_CLI_REML_REPORTER_H_
