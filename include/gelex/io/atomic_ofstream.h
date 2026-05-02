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

#ifndef GELEX_IO_ATOMIC_OFSTREAM_H_
#define GELEX_IO_ATOMIC_OFSTREAM_H_

#include <filesystem>
#include <fstream>
#include <ios>
#include <span>

namespace gelex::io::detail
{

// Writes to "<final_path>.tmp" and atomically renames to final_path on success.
// Destructor auto-commits iff stream is good() and no in-flight exception;
// otherwise removes the .tmp file.
class AtomicOfstream : public std::ofstream
{
   public:
    AtomicOfstream(
        std::filesystem::path final_path,
        std::ios::openmode mode,
        std::span<char> custom_buffer = {});

    AtomicOfstream(const AtomicOfstream&) = delete;
    AtomicOfstream(AtomicOfstream&&) = delete;
    auto operator=(const AtomicOfstream&) -> AtomicOfstream& = delete;
    auto operator=(AtomicOfstream&&) -> AtomicOfstream& = delete;

    ~AtomicOfstream() noexcept override;

    [[nodiscard]] auto final_path() const noexcept
        -> const std::filesystem::path&
    {
        return final_path_;
    }

    [[nodiscard]] auto tmp_path() const noexcept -> const std::filesystem::path&
    {
        return tmp_path_;
    }

   private:
    std::filesystem::path final_path_;
    std::filesystem::path tmp_path_;
};

}  // namespace gelex::io::detail

#endif  // GELEX_IO_ATOMIC_OFSTREAM_H_
