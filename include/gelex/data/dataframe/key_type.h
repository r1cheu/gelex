// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_DATAFRAME_KEY_TYPE_H
#define GELEX_DATA_DATAFRAME_KEY_TYPE_H

#include <cstdint>
#include <string>
#include <type_traits>

namespace gelex
{
template <typename T>
concept KeyType = std::is_arithmetic_v<T> || std::is_same_v<T, std::string>;

template <typename T>
concept ValueType
    = std::is_same_v<T, std::int32_t> || std::is_same_v<T, float>
      || std::is_same_v<T, double> || std::is_same_v<T, std::string>;

}  // namespace gelex

#endif  // GELEX_DATA_DATAFRAME_KEY_TYPE_H
