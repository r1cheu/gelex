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

#ifndef GELEX_CLI_REML_REPORTER_H_
#define GELEX_CLI_REML_REPORTER_H_

#include <memory>
#include <vector>

#include "gelex/infra/logging/reml_event.h"

namespace spdlog
{
class logger;
}

namespace gelex::cli
{

class RemlReporter
{
   public:
    RemlReporter();

    auto on_event(const RemlEmInitEvent& e) -> void;
    auto on_event(const RemlIterationEvent& e) -> void;
    auto on_event(const RemlCompleteEvent& e) const -> void;

    auto as_observer() -> RemlObserver
    {
        return [this](const RemlEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    std::shared_ptr<spdlog::logger> logger_;
    bool header_printed_ = false;
};

void print_loco_reml_summary(const std::vector<LocoRemlResult>& results);

}  // namespace gelex::cli

#endif  // GELEX_CLI_REML_REPORTER_H_
