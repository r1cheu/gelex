// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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

    return fmt::format("{}{}{}", fid, sample_id_separator, iid);
}

}  // namespace gelex

#endif  // GELEX_TESTS_SAMPLE_ID_FIXTURE_H_
