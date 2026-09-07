// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
