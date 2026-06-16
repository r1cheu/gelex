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

#ifndef APPS_CLI_REML_COMMAND_REPORTER_H_
#define APPS_CLI_REML_COMMAND_REPORTER_H_

#include <cstddef>
#include <variant>

#include "gelex/infra/logging/reml_event.h"

namespace gelex
{
class FreqModel;
class FreqState;
}  // namespace gelex

namespace CLI
{
class App;
}

namespace cli
{

class RemlCommandReporter
{
   public:
    auto show_banner() const -> void;
    auto show_config(const CLI::App& cmd) const -> void;
    auto show_result(
        const gelex::FreqModel& model,
        const gelex::FreqState& state,
        bool converged,
        size_t iter_count,
        size_t max_iter,
        double loglike) const -> void;

    auto as_observer() -> gelex::RemlObserver
    {
        return [this](const gelex::RemlEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    auto on_event(const gelex::RemlEmInitEvent& e) -> void;
    auto on_event(const gelex::RemlIterationEvent& e) -> void;

    bool header_printed_ = false;
};

}  // namespace cli

#endif  // APPS_CLI_REML_COMMAND_REPORTER_H_
