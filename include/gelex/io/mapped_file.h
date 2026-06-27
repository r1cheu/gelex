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

#ifndef GELEX_IO_MAPPED_FILE_H_
#define GELEX_IO_MAPPED_FILE_H_

#include <cstddef>
#include <memory>
#include <string>
#include <system_error>

namespace gelex::io
{

class MappedFile
{
   public:
    MappedFile() noexcept;
    ~MappedFile();
    MappedFile(MappedFile&&) noexcept;
    auto operator=(MappedFile&&) noexcept -> MappedFile&;
    MappedFile(const MappedFile&) = delete;
    auto operator=(const MappedFile&) -> MappedFile& = delete;

    auto map(const std::string& path, std::error_code& ec) -> void;

    [[nodiscard]] auto data() const noexcept -> const std::byte*;
    [[nodiscard]] auto size() const noexcept -> std::size_t;

   private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace gelex::io

#endif  // GELEX_IO_MAPPED_FILE_H_
