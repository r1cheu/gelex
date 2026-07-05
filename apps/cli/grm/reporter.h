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
#include <string_view>

#include "cli/progress_bar.h"
#include "cli/timer.h"
#include "gelex/infra/logging/grm_event.h"

namespace cli
{

class GrmReporter
{
   public:
    GrmReporter();

    static auto show_data_loaded(size_t num_samples, size_t num_snps) -> void;
    auto on_event(const gelex::GrmProgressEvent& event) -> void;
    auto finish_progress() -> void;
    static auto show_files_written(
        size_t num_files,
        std::string_view output_dir,
        std::string_view file_pattern) -> void;

    auto as_observer() -> gelex::GrmObserver
    {
        return [this](const gelex::GrmProgressEvent& e) { this->on_event(e); };
    }

   private:
    size_t progress_{0};
    cli::ProgressBar bar_;
    bool bar_active_ = false;
    gelex::SmoothEtaCalculator eta_;
    size_t global_total_ = 0;
};

}  // namespace cli

#endif  // APPS_CLI_GRM_REPORTER_H_
