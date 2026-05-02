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

#ifndef GELEX_TESTS_SAMPLE_ID_FIXTURE_H_
#define GELEX_TESTS_SAMPLE_ID_FIXTURE_H_

#include <fmt/format.h>
#include <string>
#include <string_view>

#include "gelex/data/sample_id.h"
#include "gelex/exception.h"

namespace gelex
{

inline auto make_sample_id(std::string_view fid, std::string_view iid)
    -> std::string
{
    if (fid.empty())
    {
        throw GelexException("FID cannot be empty");
    }
    if (iid.empty())
    {
        throw GelexException("IID cannot be empty");
    }

    return fmt::format("{}{}{}", fid, kSampleIdSeparator, iid);
}

}  // namespace gelex

#endif  // GELEX_TESTS_SAMPLE_ID_FIXTURE_H_
