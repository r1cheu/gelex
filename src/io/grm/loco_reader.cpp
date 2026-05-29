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

#include "gelex/io/grm/loco_reader.h"

#include <fmt/format.h>
#include <Eigen/Core>
#include <filesystem>
#include <string>

#include "gelex/data/dataframe/index.h"
#include "gelex/exception.h"
#include "gelex/io/grm/detail/reader.h"

namespace gelex
{

LocoReader::LocoReader(
    const std::filesystem::path& whole_grm_prefix,
    const dataframe::Index<std::string>& sample_index)
{
    grm::detail::GrmReader whole_reader(whole_grm_prefix);
    // Load and filter the whole GRM once during construction.
    // load_unnormalized(sample_index) returns (X_w * X_w') filtered and
    // reordered.
    g_whole_ = whole_reader.load_unnormalized(sample_index);
    // Compute trace after loading and save for LOCO calculation
    trace_whole_ = g_whole_.trace();
    k_whole_ = trace_whole_ / static_cast<double>(g_whole_.rows());
}

auto LocoReader::load_loco_grm(
    const std::filesystem::path& chr_grm_prefix,
    const dataframe::Index<std::string>& sample_index,
    Eigen::MatrixXd& target) const -> void
{
    std::filesystem::path bin_path = chr_grm_prefix.string() + ".bin";
    if (!std::filesystem::exists(bin_path))
    {
        throw GelexException(
            fmt::format(
                "LOCO error: GRM file not found: {}", bin_path.string()));
    }

    grm::detail::GrmReader chr_reader(chr_grm_prefix);

    // Load chromosome GRM filtered by the SAME sample_index to ensure
    // alignment. Use the mutable buffer to avoid reallocations.
    chr_reader.load_unnormalized(sample_index, g_chr_buffer_);

    // Compute k_i from the loaded chromosome GRM trace
    double trace_i = g_chr_buffer_.trace();
    double k_i = trace_i / static_cast<double>(g_chr_buffer_.rows());
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

    // Both g_whole_ and g_chr_buffer_ are now unnormalized (X*X') and filtered.
    target = (g_whole_ - g_chr_buffer_) / k_loco;
}

auto LocoReader::load_loco_grm(
    const std::filesystem::path& chr_grm_prefix,
    const dataframe::Index<std::string>& sample_index) const -> Eigen::MatrixXd
{
    Eigen::MatrixXd target;
    load_loco_grm(chr_grm_prefix, sample_index, target);
    return target;
}

auto LocoReader::num_samples() const noexcept -> Eigen::Index
{
    return g_whole_.rows();
}

}  // namespace gelex
