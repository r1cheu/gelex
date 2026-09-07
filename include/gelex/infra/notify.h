// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_INFRA_NOTIFY_H_
#define GELEX_INFRA_NOTIFY_H_

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
#endif  // GELEX_INFRA_NOTIFY_H_
