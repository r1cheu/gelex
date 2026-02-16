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

#ifndef GELEX_INTERNAL_DATA_GRM_ID_WRITER_H_
#define GELEX_INTERNAL_DATA_GRM_ID_WRITER_H_

#include <filesystem>
#include <fstream>
#include <span>
#include <string>
#include <string_view>
#include <utility>

namespace gelex::detail
{

class GrmIdWriter
{
   public:
    explicit GrmIdWriter(const std::filesystem::path& file_path);

    GrmIdWriter(const GrmIdWriter&) = delete;
    GrmIdWriter(GrmIdWriter&&) noexcept = default;
    GrmIdWriter& operator=(const GrmIdWriter&) = delete;
    GrmIdWriter& operator=(GrmIdWriter&&) noexcept = default;
    ~GrmIdWriter() = default;

    auto write(std::span<const std::string> ids) -> void;

    [[nodiscard]] auto path() const noexcept -> const std::filesystem::path&
    {
        return path_;
    }

   private:
    std::filesystem::path path_;
    std::ofstream file_;

    static auto split_id(std::string_view id)
        -> std::pair<std::string_view, std::string_view>;
};

}  // namespace gelex::detail

#endif  // GELEX_INTERNAL_DATA_GRM_ID_WRITER_H_
