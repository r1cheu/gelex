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

#ifndef GELEX_DATA_DATAFRAME_CONSTANTS_H
#define GELEX_DATA_DATAFRAME_CONSTANTS_H

#include <array>
#include <string_view>

namespace gelex::df
{

inline constexpr char kSeparator = '\x1F';

inline constexpr std::array kDefaultNaRep = {
    std::string_view{"NA"},
    std::string_view{"NaN"},
    std::string_view{"nan"},
    std::string_view{"null"},
    std::string_view{"NULL"},
    std::string_view{"."},
};

}  // namespace gelex::df

#endif  // GELEX_DATA_DATAFRAME_CONSTANTS_H
