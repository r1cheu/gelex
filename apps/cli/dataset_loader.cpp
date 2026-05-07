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

#include "dataset_loader.h"

#include <fmt/format.h>
#include <utility>

#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"

namespace gelex::cli
{

auto intersect_or_throw(
    std::vector<const dataframe::Index<std::string>*> indices,
    const DatasetObserver& observer,
    std::string_view what) -> dataframe::Index<std::string>
{
    auto common = dataframe::intersect<std::string>(indices);
    notify(observer, IntersectionEvent{.common_samples = common.size()});

    if (common.size() == 0)
    {
        throw GelexException(
            fmt::format(
                "No common samples across {}. Check that sample IDs match "
                "across input files.",
                what));
    }
    return common;
}

}  // namespace gelex::cli
