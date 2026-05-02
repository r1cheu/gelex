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

#ifndef GELEX_IO_LOCISTATS_WRITER_H_
#define GELEX_IO_LOCISTATS_WRITER_H_

#include <cstdint>
#include <span>
#include <string_view>

#include <Eigen/Core>

#include "gelex/io/detail/binary_writer.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class LociStatsWriter
{
   public:
    explicit LociStatsWriter(std::string_view output_path);

    auto write(
        EffectType effect,
        uint8_t method,
        const Eigen::VectorXd& mean,
        const Eigen::VectorXd* stddev = nullptr,
        std::span<const int64_t> mono_indices = {}) -> void;

   private:
    io::detail::BinaryWriter writer_;
};

}  // namespace gelex

#endif  // GELEX_IO_LOCISTATS_WRITER_H_
