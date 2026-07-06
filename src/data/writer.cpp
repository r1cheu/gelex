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

#include "gelex/data/writer.h"

#include <ios>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/data/sample_id.h"
#include "gelex/exception.h"
#include "gelex/io/detail/atomic_output_stream.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex
{
auto write_grm_ids(const std::string& prefix, std::span<const std::string> ids)
    -> void
{
    detail::TextWriter writer(prefix + ".id");
    for (const auto& id : ids)
    {
        auto [fid, iid] = split_sample_id(id);
        writer.write(fmt::format("{}\t{}", fid, iid));
    }
}

auto write_grm(
    const std::string& prefix,
    const Eigen::Ref<const Eigen::MatrixXd>& grm,
    std::span<const std::string> ids) -> void
{
    detail::AtomicOutputStream file(prefix + ".bin", std::ios::binary);

    auto n = grm.rows();
    auto m = grm.cols();

    if (n != m)
    {
        throw GelexException(
            fmt::format(
                "{}: GRM must be square, got {}x{}",
                file.path().string(),
                n,
                m));
    }
    if (n != static_cast<Eigen::Index>(ids.size()))
    {
        throw GelexException(
            fmt::format(
                "{}: Number of IDs ({}) does not match GRM size ({})",
                file.path().string(),
                ids.size(),
                n));
    }

    std::vector<float> row_buffer(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i)
    {
        for (Eigen::Index j = 0; j <= i; ++j)
        {
            row_buffer[static_cast<size_t>(j)] = static_cast<float>(grm(i, j));
        }
        file.write(
            reinterpret_cast<const char*>(row_buffer.data()),
            static_cast<std::streamsize>((i + 1) * sizeof(float)));
    }
    file.commit();
    write_grm_ids(prefix, ids);
}

}  // namespace gelex
