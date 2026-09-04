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

#ifndef APPS_CLI_GRM_REPORTER_H_
#define APPS_CLI_GRM_REPORTER_H_

#include <cstddef>

#include "gelex/data/grm/progress.h"

#include "cli/progress.h"

namespace cli
{

class GrmReporter
{
   public:
    explicit GrmReporter(std::size_t total);

    static auto show_data_loaded(size_t num_samples, size_t num_snps) -> void;
    auto on_event(const gelex::GrmProgressEvent& event) -> void;

    auto as_observer() -> gelex::GrmObserver
    {
        return [this](const gelex::GrmProgressEvent& e) { this->on_event(e); };
    }

   private:
    cli::Progress progress_;
    decltype(cli::make_rate()) estimate_rate_;
    decltype(cli::make_eta(std::size_t{})) estimate_eta_;
};

}  // namespace cli

#endif  // APPS_CLI_GRM_REPORTER_H_
