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

#include "gelex/io/locistats_reader.h"

#include <cstdint>
#include <string_view>

namespace gelex
{

LociStatsReader::LociStatsReader(std::string_view file_path)
    : reader_(file_path)
{
}

auto LociStatsReader::has(EffectType effect) const -> bool
{
    return reader_.has_section(effect, detail::DataKind::LociStats);
}

auto LociStatsReader::read(EffectType effect) const -> LociStats
{
    auto stats_map = reader_.map<double>(effect, detail::DataKind::LociStats);

    LociStats data;
    data.mean = stats_map.col(0);

    if (reader_.has_section(effect, detail::DataKind::GenoMethod))
    {
        auto method_map
            = reader_.map<uint8_t>(effect, detail::DataKind::GenoMethod);
        data.method = static_cast<GenotypeProcessMethod>(method_map(0, 0));
    }

    if (stats_map.cols() == 2)
    {
        data.stddev = Eigen::VectorXd(stats_map.col(1));
    }

    if (reader_.has_section(effect, detail::DataKind::MonoIndices))
    {
        auto mono_mat
            = reader_.map<int64_t>(effect, detail::DataKind::MonoIndices);
        const auto* src = mono_mat.col(0).data();
        data.mono_indices.assign(src, src + mono_mat.rows());
    }

    return data;
}

}  // namespace gelex
