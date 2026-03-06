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

#ifndef GELEX_INFRA_LOGGING_NOTIFY_H_
#define GELEX_INFRA_LOGGING_NOTIFY_H_

#include <utility>

namespace gelex
{

template <typename Observer, typename Event>
auto notify(const Observer& observer, Event&& event) -> void
{
    if (observer)
    {
        observer(std::forward<Event>(event));
    }
}

}  // namespace gelex
#endif  // GELEX_INFRA_LOGGING_NOTIFY_H_
