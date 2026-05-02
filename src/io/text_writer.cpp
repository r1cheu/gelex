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

namespace gelex::io::detail
{

TextWriter::TextWriter(const std::filesystem::path& path)
    : buf_{}, ofs_(path, std::ios::out, buf_)
{
}

auto TextWriter::write_header(std::initializer_list<std::string_view> columns)
    -> void
{
    bool first = true;
    for (const auto col : columns)
    {
        if (!first)
        {
            ofs_.put('\t');
        }
        ofs_.write(col.data(), static_cast<std::streamsize>(col.size()));
        first = false;
    }
    ofs_.put('\n');
}

auto TextWriter::write(std::string_view line) -> void
{
    ofs_.write(line.data(), static_cast<std::streamsize>(line.size()));
    ofs_.put('\n');
}

auto TextWriter::path() const noexcept -> const std::filesystem::path&
{
    return ofs_.final_path();
}

}  // namespace gelex::io::detail
