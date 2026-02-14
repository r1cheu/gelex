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

#ifndef GELEX_INTERNAL_DATAFRAME_TEXT_UTILS_H_
#define GELEX_INTERNAL_DATAFRAME_TEXT_UTILS_H_

#include <string_view>
#include <vector>

namespace gelex::detail
{

auto split_line_preserve_empty(std::string_view line, char delimiter)
    -> std::vector<std::string_view>;

auto detect_delimiter(std::string_view probe) -> char;

}  // namespace gelex::detail

#endif  // GELEX_INTERNAL_DATAFRAME_TEXT_UTILS_H_
