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

#ifndef GELEX_DATA_DATAFRAME_KEY_TYPE_H
#define GELEX_DATA_DATAFRAME_KEY_TYPE_H

#include <cstdint>
#include <string>
#include <type_traits>

namespace gelex::df
{
template <typename T>
concept KeyType = std::is_arithmetic_v<T> || std::is_same_v<T, std::string>;

template <typename T>
concept ValueType
    = std::is_same_v<T, std::int32_t> || std::is_same_v<T, float>
      || std::is_same_v<T, double> || std::is_same_v<T, std::string>;

}  // namespace gelex::df

#endif  // GELEX_DATA_DATAFRAME_KEY_TYPE_H
