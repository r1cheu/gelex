// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/io/detail/parser.h"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <vector>

namespace gelex::detail
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

}  // namespace gelex::detail
