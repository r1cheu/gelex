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

#include "gelex/io/sbin_writer.h"

#include <cstdint>
#include <format>
#include <span>
#include <string_view>

#include "gelex/exception.h"

namespace gelex
{

SbinWriter::SbinWriter(std::string_view output_path) : writer_(output_path) {}

auto SbinWriter::write(
    detail::EffectType effect,
    const Eigen::VectorXd& mean,
    const Eigen::VectorXd* stddev,
    std::span<const int64_t> mono_indices) -> void
{
    if (stddev != nullptr && mean.size() != stddev->size())
    {
        throw ArgumentValidationException(
            std::format(
                "SbinWriter: mean size ({}) != stddev size ({})",
                mean.size(),
                stddev->size()));
    }

    const auto n_snps = static_cast<uint64_t>(mean.size());
    const uint64_t n_cols = (stddev != nullptr) ? 2 : 1;

    auto stats_handle = writer_.reserve(
        {effect, detail::DataKind::SnpStats},
        detail::binary_format::dtype_code<double>(),
        n_snps,
        n_cols);
    writer_.write(
        stats_handle,
        reinterpret_cast<const char*>(mean.data()),
        static_cast<std::streamsize>(mean.size() * sizeof(double)));

    if (stddev != nullptr)
    {
        writer_.write(
            stats_handle,
            reinterpret_cast<const char*>(stddev->data()),
            static_cast<std::streamsize>(stddev->size() * sizeof(double)));
    }

    if (!mono_indices.empty())
    {
        auto mono_handle = writer_.reserve(
            {effect, detail::DataKind::MonoIndices},
            detail::binary_format::dtype_code<int64_t>(),
            static_cast<uint64_t>(mono_indices.size()),
            1);
        writer_.write(
            mono_handle,
            reinterpret_cast<const char*>(mono_indices.data()),
            static_cast<std::streamsize>(
                mono_indices.size() * sizeof(int64_t)));
    }
}

auto SbinWriter::finalize() -> void
{
    writer_.finalize();
}

}  // namespace gelex
