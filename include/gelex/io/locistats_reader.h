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

#ifndef GELEX_IO_LOCISTATS_READER_H_
#define GELEX_IO_LOCISTATS_READER_H_

#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

#include <Eigen/Core>

#include "gelex/io/detail/binary_reader.h"
#include "gelex/types/genetic_effect_type.h"
#include "gelex/types/genotype_process_method.h"

namespace gelex
{

struct LociStats
{
    GenotypeProcessMethod method{};
    Eigen::VectorXd mean;
    std::optional<Eigen::VectorXd> stddev;
    std::vector<int64_t> mono_indices;
};

class LociStatsReader
{
   public:
    explicit LociStatsReader(std::string_view file_path);

    [[nodiscard]] auto read(EffectType effect) const -> LociStats;
    [[nodiscard]] auto has(EffectType effect) const -> bool;

   private:
    io::detail::BinaryReader reader_;
};

}  // namespace gelex

#endif  // GELEX_IO_LOCISTATS_READER_H_
