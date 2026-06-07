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

#include "gelex/io/detail/text_writer.h"
#include <exception>
#include <filesystem>
#include <initializer_list>
#include <ios>
#include <string_view>

namespace gelex::io::detail
{

TextWriter::TextWriter(const std::filesystem::path& path)
    : ofs_(path, std::ios::out)
{
}

TextWriter::~TextWriter() noexcept
{
    if (std::uncaught_exceptions() > 0)
    {
        return;
    }
    try
    {
        ofs_.commit();
    }
    catch (...)  // NOLINT(bugprone-empty-catch): dtor must be noexcept
    {
    }
}

auto TextWriter::write_header(std::initializer_list<std::string_view> columns)
    -> void
{
    bool first = true;
    for (const auto col : columns)
    {
        if (!first)
        {
            ofs_.write("\t");
        }
        ofs_.write(col);
        first = false;
    }
    ofs_.write("\n");
}

auto TextWriter::write(std::string_view line) -> void
{
    ofs_.write(line);
    ofs_.write("\n");
}

auto TextWriter::path() const noexcept -> const std::filesystem::path&
{
    return ofs_.path();
}

}  // namespace gelex::io::detail
