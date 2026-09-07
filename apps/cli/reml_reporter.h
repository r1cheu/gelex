// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_REML_REPORTER_H_
#define APPS_CLI_REML_REPORTER_H_

#include <cstddef>
#include <variant>
#include <vector>

#include "gelex/freq/reml/progress.h"
#include "gelex/freq/reml/summary.h"

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
    bool has_prev_ = false;
    double prev_loglike_ = 0.0;
    std::vector<double> prev_variances_;
};

void print_loco_reml_summary(const std::vector<gelex::LocoRemlResult>& results);

}  // namespace cli

#endif  // APPS_CLI_REML_REPORTER_H_
