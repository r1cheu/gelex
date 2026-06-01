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

#ifndef GELEX_IO_GRM_WRITER_H_
#define GELEX_IO_GRM_WRITER_H_

#include <cstddef>
#include <filesystem>
#include <span>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/io/detail/atomic_ofstream.h"

namespace gelex
{

class GrmBinWriter
{
   public:
    static constexpr size_t DEFAULT_BUFFER_SIZE
        = static_cast<size_t>(64 * 1024);

    explicit GrmBinWriter(const std::filesystem::path& file_path);

    GrmBinWriter(const GrmBinWriter&) = delete;
    GrmBinWriter(GrmBinWriter&&) = delete;
    auto operator=(const GrmBinWriter&) -> GrmBinWriter& = delete;
    auto operator=(GrmBinWriter&&) -> GrmBinWriter& = delete;
    ~GrmBinWriter() = default;

    auto write(const Eigen::Ref<const Eigen::MatrixXd>& grm) -> void;

    [[nodiscard]] auto path() const noexcept -> const std::filesystem::path&
    {
        return file_.final_path();
    }

   private:
    std::vector<char> io_buffer_;
    io::detail::AtomicOfstream file_;
};

auto write_grm_ids(
    const std::filesystem::path& file_path,
    std::span<const std::string> ids) -> void;

}  // namespace gelex

#endif  // GELEX_IO_GRM_WRITER_H_
