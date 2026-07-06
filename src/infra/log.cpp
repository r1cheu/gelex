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

#include "gelex/infra/log.h"

#include <string_view>
#include <utility>

namespace gelex
{

namespace
{

auto sink() -> Sink&
{
    static Sink instance;
    return instance;
}

}  // namespace

void set_sink(Sink s)
{
    sink() = std::move(s);
}

void info(std::string_view message)
{
    if (auto& s = sink())
    {
        s(Level::Info, message);
    }
}

void warn(std::string_view message)
{
    if (auto& s = sink())
    {
        s(Level::Warn, message);
    }
}

void error(std::string_view message)
{
    if (auto& s = sink())
    {
        s(Level::Error, message);
    }
}

}  // namespace gelex
