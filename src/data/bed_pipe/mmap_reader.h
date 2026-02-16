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

#ifndef GELEX_DATA_BED_PIPE_MMAP_READER_H_
#define GELEX_DATA_BED_PIPE_MMAP_READER_H_

#include <cstddef>
#include <filesystem>

#include <Eigen/Core>

#include "mio.h"

namespace gelex::detail
{

class BedMmapReader
{
   public:
    BedMmapReader(
        const std::filesystem::path& bed_path,
        Eigen::Index num_raw_snps,
        Eigen::Index bytes_per_variant);

    [[nodiscard]] auto chunk_ptr(Eigen::Index start_snp, Eigen::Index num_snps)
        const -> const uint8_t*;

    [[nodiscard]] auto bytes_per_variant() const -> Eigen::Index
    {
        return bytes_per_variant_;
    }

   private:
    mio::mmap_source mmap_;
    const uint8_t* payload_ptr_ = nullptr;
    size_t payload_num_snps_ = 0;
    Eigen::Index bytes_per_variant_ = 0;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_BED_PIPE_MMAP_READER_H_
