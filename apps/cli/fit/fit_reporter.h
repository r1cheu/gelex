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

#ifndef GELEX_CLI_FIT_REPORTER_H_
#define GELEX_CLI_FIT_REPORTER_H_

#include <atomic>
#include <memory>
#include <string>

#include "gelex/infra/detail/indicator.h"

namespace gelex
{
struct FitConfigLoadedEvent;
struct FitMcmcProgressEvent;
struct FitResultsSavedEvent;
}  // namespace gelex

namespace spdlog
{
class logger;
}

namespace gelex::cli
{

class FitReporter
{
   public:
    FitReporter();

    auto on_event(const FitConfigLoadedEvent& event) const -> void;
    auto on_event(const FitMcmcProgressEvent& event) -> void;
    auto on_event(const FitResultsSavedEvent& event) const -> void;

   private:
    std::shared_ptr<spdlog::logger> logger_;
    std::atomic<size_t> iter_{0};
    detail::ProgressBar bar_;
    bool init_progress_ = false;
    std::string stats_;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_FIT_REPORTER_H_
