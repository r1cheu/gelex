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

#ifndef GELEX_IO_BINARY_WRITER_H_
#define GELEX_IO_BINARY_WRITER_H_

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string_view>
#include <vector>

#include "gelex/io/binary_format.h"

namespace gelex::detail
{

class BinaryWriter
{
   public:
    struct ReservedSection
    {
        TocEntry entry;
        uint64_t cursor{0};
    };

    explicit BinaryWriter(std::string_view output_path);

    BinaryWriter(const BinaryWriter&) = delete;
    BinaryWriter(BinaryWriter&&) noexcept = default;
    auto operator=(const BinaryWriter&) -> BinaryWriter& = delete;
    auto operator=(BinaryWriter&&) noexcept -> BinaryWriter& = default;
    ~BinaryWriter() = default;

    auto reserve(SectionKey key, uint8_t dtype, uint64_t rows, uint64_t cols)
        -> size_t;
    auto write(size_t handle, const char* data, std::streamsize bytes) -> void;
    auto finalize() -> void;

   private:
    auto open() -> std::ofstream&;
    auto check_duplicate_key(const SectionKey& key) -> void;

    std::filesystem::path output_path_;
    std::vector<ReservedSection> reserved_;

    std::ofstream file_;
    uint64_t next_offset_{0};
};

}  // namespace gelex::detail

#endif  // GELEX_IO_BINARY_WRITER_H_
