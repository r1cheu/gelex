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

#ifndef GELEX_IO_TEXT_WRITER_H_
#define GELEX_IO_TEXT_WRITER_H_

#include <array>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <string_view>

namespace gelex::detail
{

class TextWriter
{
   public:
    static constexpr std::streamsize kBufSize = 64 * 1024;

    explicit TextWriter(const std::filesystem::path& path);
    TextWriter(const TextWriter&) = delete;
    TextWriter(TextWriter&&) = delete;
    auto operator=(const TextWriter&) -> TextWriter& = delete;
    auto operator=(TextWriter&&) -> TextWriter& = delete;
    ~TextWriter() noexcept;

    auto write_header(std::initializer_list<std::string_view> columns) -> void;
    auto write(std::string_view line) -> void;

    [[nodiscard]] auto path() const noexcept -> const std::filesystem::path&;

   private:
    std::filesystem::path path_;
    std::array<char, kBufSize> buf_;
    std::ofstream ofs_;
};

}  // namespace gelex::detail

#endif  // GELEX_IO_TEXT_WRITER_H_
