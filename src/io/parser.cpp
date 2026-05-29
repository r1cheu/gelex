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

#include "gelex/io/detail/parser.h"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <vector>

namespace gelex::io::detail
{

size_t count_total_lines(const std::filesystem::path& path)
{
    constexpr auto buffer_size = static_cast<size_t>(1024 * 128);

    std::vector<char> buffer(buffer_size);

    auto file = open_file<std::ifstream>(path, std::ios::in | std::ios::binary);

    file.rdbuf()->pubsetbuf(buffer.data(), buffer_size);

    size_t line_count = 0;

    while (file)
    {
        file.read(buffer.data(), buffer_size);
        std::streamsize count = file.gcount();
        if (count == 0)
        {
            break;
        }

        line_count
            += std::ranges::count(buffer.begin(), buffer.begin() + count, '\n');

        if (file.eof() && count > 0 && buffer[count - 1] != '\n')
        {
            line_count++;
        }
    }

    return line_count;
}

}  // namespace gelex::io::detail
