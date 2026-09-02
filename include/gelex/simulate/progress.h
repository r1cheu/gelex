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
#ifndef GELEX_SIMULATE_PROGRESS_H_
#define GELEX_SIMULATE_PROGRESS_H_

#include <cstddef>
#include <functional>

namespace gelex
{

struct SimulateProgressEvent
{
    std::size_t total;
    std::size_t current;
    bool done;
};

using SimulateObserver
    = std::function<void(const SimulateProgressEvent& event)>;

}  // namespace gelex

#endif  // GELEX_SIMULATE_PROGRESS_H_
