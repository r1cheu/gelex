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

#ifndef GELEX_IO_LOCO_READER_H_
#define GELEX_IO_LOCO_READER_H_

#include <filesystem>
#include <string>

#include <Eigen/Core>

#include "gelex/data/dataframe/index.h"

namespace gelex
{

class LocoReader
{
   public:
    LocoReader(
        const std::filesystem::path& whole_grm_prefix,
        const dataframe::Index<std::string>& sample_index);

    auto load_into(
        const std::filesystem::path& chr_grm_prefix,
        const dataframe::Index<std::string>& sample_index,
        Eigen::MatrixXd& target) const -> void;

   private:
    Eigen::MatrixXd g_whole_;
    double k_whole_{};
    double trace_whole_{};
};

}  // namespace gelex

#endif  // GELEX_IO_LOCO_READER_H_
