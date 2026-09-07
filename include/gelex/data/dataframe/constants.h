// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_DATAFRAME_CONSTANTS_H
#define GELEX_DATA_DATAFRAME_CONSTANTS_H

#include <array>
#include <string_view>

namespace gelex
{

inline constexpr char separator = '\x1F';

inline constexpr std::string_view intercept_name = "Intercept";

inline constexpr std::array default_na_rep = {
    std::string_view{""},
    std::string_view{"NA"},
    std::string_view{"NaN"},
    std::string_view{"nan"},
    std::string_view{"null"},
    std::string_view{"NULL"},
    std::string_view{"."},
};

}  // namespace gelex

#endif  // GELEX_DATA_DATAFRAME_CONSTANTS_H
