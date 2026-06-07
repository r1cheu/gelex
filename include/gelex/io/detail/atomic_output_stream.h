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

#ifndef GELEX_IO_DETAIL_ATOMIC_OUTPUT_STREAM_H_
#define GELEX_IO_DETAIL_ATOMIC_OUTPUT_STREAM_H_

#include <filesystem>
#include <fstream>
#include <ios>
#include <string_view>

namespace gelex::io::detail
{

class AtomicOutputStream
{
   public:
    AtomicOutputStream(std::filesystem::path path, std::ios::openmode mode);

    AtomicOutputStream(const AtomicOutputStream&) = delete;
    AtomicOutputStream(AtomicOutputStream&&) = delete;
    auto operator=(const AtomicOutputStream&) -> AtomicOutputStream& = delete;
    auto operator=(AtomicOutputStream&&) -> AtomicOutputStream& = delete;

    ~AtomicOutputStream() noexcept;

    auto write(const char* data, std::streamsize size) -> void;
    auto write(std::string_view text) -> void;
    auto seek(std::streamoff offset) -> void;
    auto commit() -> void;

    [[nodiscard]] auto path() const noexcept -> const std::filesystem::path&
    {
        return path_;
    }

   private:
    std::filesystem::path path_;
    std::filesystem::path tmp_path_;
    std::ofstream file_;
    bool committed_{false};
};

}  // namespace gelex::io::detail

#endif  // GELEX_IO_DETAIL_ATOMIC_OUTPUT_STREAM_H_
