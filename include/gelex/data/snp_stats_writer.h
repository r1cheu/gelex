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

#ifndef GELEX_DATA_SNP_STATS_WRITER_H_
#define GELEX_DATA_SNP_STATS_WRITER_H_

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <string_view>
#include <vector>

namespace gelex::detail
{

class SnpStatsWriter
{
   public:
    static constexpr auto kDefaultBufferSize = static_cast<size_t>(64 * 1024);

    explicit SnpStatsWriter(const std::filesystem::path& file_path);

    SnpStatsWriter(const SnpStatsWriter&) = delete;
    SnpStatsWriter(SnpStatsWriter&&) noexcept = default;
    SnpStatsWriter& operator=(const SnpStatsWriter&) = delete;
    SnpStatsWriter& operator=(SnpStatsWriter&&) noexcept = default;
    ~SnpStatsWriter() = default;

    void write(
        int64_t num_samples,
        std::span<const int64_t> monomorphic_indices,
        std::span<const double> means,
        std::span<const double> stddevs);

    [[nodiscard]] auto path() const noexcept -> const std::filesystem::path&
    {
        return path_;
    }

   private:
    void write_data(
        const void* data,
        std::streamsize size,
        std::string_view error_msg);

    template <typename T>
    void write_data(std::span<const T> data, std::string_view error_msg)
    {
        write_data(
            data.data(),
            static_cast<std::streamsize>(data.size_bytes()),
            error_msg);
    }

    static void check_monomorphic_indices(
        std::span<const int64_t> monomorphic_indices,
        int64_t num_variants);

    std::filesystem::path path_;
    std::vector<char> io_buffer_;
    std::ofstream file_;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_SNP_STATS_WRITER_H_
