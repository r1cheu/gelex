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

#ifndef GELEX_INFRA_LOG_H_
#define GELEX_INFRA_LOG_H_

#include <cstdint>
#include <functional>
#include <string_view>

namespace gelex
{

enum class Level : std::uint8_t
{
    Info,
    Warn,
    Error
};

using Sink = std::function<void(Level, std::string_view)>;

void set_sink(Sink sink);

void info(std::string_view message);
void warn(std::string_view message);
void error(std::string_view message);

}  // namespace gelex

#endif  // GELEX_INFRA_LOG_H_
