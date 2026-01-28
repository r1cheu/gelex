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

#include "binary_matrix_writer.h"

#include "gelex/exception.h"
#include "parser.h"

namespace gelex::detail
{

BinaryMatrixWriter::BinaryMatrixWriter(const std::filesystem::path& file_path)
    : path_(file_path), io_buffer_(kDefaultBufferSize)
{
    file_ = detail::open_file<std::ofstream>(
        path_, std::ios::binary | std::ios::trunc, io_buffer_);
}

void BinaryMatrixWriter::write(const Eigen::Ref<const Eigen::MatrixXd>& matrix)
{
    if (matrix.size() == 0)
    {
        return;
    }

    file_.write(
        reinterpret_cast<const char*>(matrix.data()),
        static_cast<std::streamsize>(matrix.size() * sizeof(double)));

    if (!file_.good())
    {
        throw FileWriteException(
            std::format(
                "{}: failed to write matrix data to binary file",
                path_.string()));
    }
}

}  // namespace gelex::detail
