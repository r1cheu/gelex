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

#include "gelex/io/grm/writer.h"

#include <fmt/format.h>
#include <Eigen/Core>
#include <filesystem>
#include <ios>
#include <span>
#include <string>

#include "gelex/data/sample_id.h"
#include "gelex/exception.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex
{

GrmBinWriter::GrmBinWriter(const std::filesystem::path& file_path)
    : file_(file_path, std::ios::binary | std::ios::trunc)
{
}

auto GrmBinWriter::write(const Eigen::Ref<const Eigen::MatrixXd>& grm) -> void
{
    const Eigen::Index n = grm.rows();
    if (n == 0)
    {
        file_.commit();
        return;
    }

    if (grm.rows() != grm.cols())
    {
        throw GelexException(
            fmt::format(
                "{}: GRM must be square, got {}x{}",
                file_.path().string(),
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

    file_.commit();
}

auto write_grm_ids(
    const std::filesystem::path& file_path,
    std::span<const std::string> ids) -> void
{
    io::detail::TextWriter writer(file_path);
    for (const auto& id : ids)
    {
        auto [fid, iid] = split_sample_id(id);
        writer.write(fmt::format("{}\t{}", fid, iid));
    }
}

}  // namespace gelex
