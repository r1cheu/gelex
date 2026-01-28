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

#include "grm_bin_writer.h"

#include <format>

#include "gelex/exception.h"
#include "parser.h"

namespace gelex::detail
{

GrmBinWriter::GrmBinWriter(const std::filesystem::path& file_path)
    : path_(file_path), io_buffer_(kDefaultBufferSize)
{
    file_ = detail::open_file<std::ofstream>(
        path_, std::ios::binary | std::ios::trunc, io_buffer_);
}

auto GrmBinWriter::write(const Eigen::Ref<const Eigen::MatrixXd>& grm) -> void
{
    const Eigen::Index n = grm.rows();
    if (n == 0)
    {
        return;
    }

    if (grm.rows() != grm.cols())
    {
        throw InvalidInputException(
            std::format(
                "{}: GRM must be square, got {}x{}",
                path_.string(),
                n,
                grm.cols()));
    }

    // Write lower triangle (including diagonal) as float32
    // Order: (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), ...
    for (Eigen::Index i = 0; i < n; ++i)
    {
        for (Eigen::Index j = 0; j <= i; ++j)
        {
            auto value = static_cast<float>(grm(i, j));
            file_.write(reinterpret_cast<const char*>(&value), sizeof(float));
        }
    }

    if (!file_.good())
    {
        throw FileWriteException(
            std::format(
                "{}: failed to write GRM data to binary file", path_.string()));
    }
}

}  // namespace gelex::detail
