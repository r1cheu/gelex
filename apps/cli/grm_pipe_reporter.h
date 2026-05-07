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

#ifndef GELEX_CLI_GRM_PIPE_REPORTER_H_
#define GELEX_CLI_GRM_PIPE_REPORTER_H_

#include <memory>
#include <variant>

#include "gelex/infra/logging/grm_pipe_event.h"

namespace spdlog
{
class logger;
}

namespace gelex::cli
{

class GrmPipeReporter
{
   public:
    GrmPipeReporter();

    auto on_event(const GrmLoadedEvent& event) const -> void;

    auto as_observer() -> GrmPipeObserver
    {
        return [this](const GrmPipeEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    std::shared_ptr<spdlog::logger> logger_;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_GRM_PIPE_REPORTER_H_
