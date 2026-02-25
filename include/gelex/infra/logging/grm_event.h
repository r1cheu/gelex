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

#ifndef GELEX_INFRA_LOGGING_GRM_EVENT_H_
#define GELEX_INFRA_LOGGING_GRM_EVENT_H_

#include <functional>
#include <string>
#include <variant>
#include <vector>

#include "gelex/types/freq_effect.h"

namespace gelex
{

struct GrmConfigLoadedEvent
{
    std::string method;
    gelex::freq::GrmType mode;
    bool do_loco;
    int chunk_size;
    int threads;
};

struct GrmDataLoadedEvent
{
    size_t num_samples;
    size_t num_snps;
};

struct GrmProgressEvent
{
    size_t current;
    size_t total;
    bool done;
};

struct GrmFilesWrittenEvent
{
    std::vector<std::string> file_paths;
    std::string time_elapsed;
    std::string output_dir;
    std::string file_pattern;
};

using GrmEvent = std::variant<
    GrmConfigLoadedEvent,
    GrmDataLoadedEvent,
    GrmProgressEvent,
    GrmFilesWrittenEvent>;

using GrmObserver = std::function<void(const GrmEvent&)>;

}  // namespace gelex

#endif  // GELEX_INFRA_LOGGING_GRM_EVENT_H_
