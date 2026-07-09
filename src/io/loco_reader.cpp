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

#include "gelex/io/loco_reader.h"

#include <fmt/format.h>
#include <Eigen/Core>
#include <filesystem>
#include <string>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"

namespace gelex
{

LocoReader::LocoReader(const Eigen::MatrixXd& whole_grm)
    : g_whole_(&whole_grm),
      k_whole_(whole_grm.trace() / static_cast<double>(whole_grm.rows()))
{
}

auto LocoReader::load_into(
    const std::filesystem::path& chr_grm_prefix,
    const DataFrameIndex<std::string>& sample_index,
    Eigen::MatrixXd& target) const -> void
{
    std::filesystem::path bin_path = chr_grm_prefix.string() + ".bin";
    if (!std::filesystem::exists(bin_path))
    {
        throw GelexException(
            fmt::format(
                "LOCO error: GRM file not found: {}", bin_path.string()));
    }

    Eigen::MatrixXd g_chr
        = read_grm(chr_grm_prefix.string(), &sample_index, false);

    double trace_i = g_chr.trace();
    double k_i = trace_i / static_cast<double>(g_chr.rows());
    double k_loco = k_whole_ - k_i;
    if (k_loco <= 0)
    {
        throw GelexException(
            fmt::format(
                "LOCO error: Chromosome GRM denominator ({}) is greater than "
                "or equal to Whole GRM denominator ({})",
                k_i,
                k_whole_));
    }

    target = (*g_whole_ - g_chr) / k_loco;
}

}  // namespace gelex
