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

#ifndef GELEX_DATA_GRM_BIN_WRITER_H
#define GELEX_DATA_GRM_BIN_WRITER_H

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <vector>

#include <Eigen/Core>

namespace gelex::detail
{

class GrmBinWriter
{
   public:
    static constexpr size_t kDefaultBufferSize = static_cast<size_t>(64 * 1024);

    explicit GrmBinWriter(const std::filesystem::path& file_path);

    GrmBinWriter(const GrmBinWriter&) = delete;
    GrmBinWriter(GrmBinWriter&&) noexcept = default;
    GrmBinWriter& operator=(const GrmBinWriter&) = delete;
    GrmBinWriter& operator=(GrmBinWriter&&) noexcept = default;
    ~GrmBinWriter() = default;

    // Write GRM matrix (unnormalized)
    // Format: [float32 lower triangle]
    // Order: (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), ...
    auto write(const Eigen::Ref<const Eigen::MatrixXd>& grm) -> void;

    [[nodiscard]] auto path() const noexcept -> const std::filesystem::path&
    {
        return path_;
    }

   private:
    std::filesystem::path path_;
    std::vector<char> io_buffer_;
    std::ofstream file_;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_GRM_BIN_WRITER_H
