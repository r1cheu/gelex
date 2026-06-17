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

#include "gelex/io/gwas/joint_cov_writer.h"

#include <cstddef>
#include <cstdint>
#include <exception>
#include <ios>
#include <iterator>
#include <string>
#include <string_view>

#include <fmt/base.h>
#include <fmt/compile.h>

#include "gelex/algo/gwas/assoc_tester.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/exception.h"

namespace gelex::gwas
{

constexpr size_t JOINT_COV_BUFFER_FLUSH_THRESHOLD
    = static_cast<size_t>(64 * 1024);

JointCovWriter::JointCovWriter(
    std::string_view out_prefix,
    const dataframe::DataFrame<std::string>& bim)
    : keys_(bim.index().keys()),
      chrom_(bim["chrom"].as<std::string>()),
      pos_(bim["pos"].as<std::int32_t>()),
      a1_(bim["A1"].as<std::string>()),
      a2_(bim["A2"].as<std::string>()),
      ofs_(std::string(out_prefix) + ".cov", std::ios::out | std::ios::binary)
{
    line_buffer_.reserve(JOINT_COV_BUFFER_FLUSH_THRESHOLD);

    fmt::format_to(
        std::back_inserter(line_buffer_),
        FMT_COMPILE(
            "CHR\tSNP\tBP\tA1\tA2\tA1FREQ\tBETA_A\tBETA_D\t"
            "VAR_A\tVAR_D\tCOV_AD\n"));

    ofs_.write(
        line_buffer_.data(), static_cast<std::streamsize>(line_buffer_.size()));
    line_buffer_.clear();
}

JointCovWriter::~JointCovWriter() noexcept
{
    if (std::uncaught_exceptions() > 0)
    {
        return;
    }
    try
    {
        if (line_buffer_.size() != 0U)
        {
            ofs_.write(
                line_buffer_.data(),
                static_cast<std::streamsize>(line_buffer_.size()));
            line_buffer_.clear();
        }
        ofs_.commit();
    }
    catch (...)  // NOLINT(bugprone-empty-catch): dtor must be noexcept
    {
    }
}

auto JointCovWriter::write(std::size_t start, const TestResults& results)
    -> void
{
    if (!results.dominance || !results.joint_covariance)
    {
        throw GelexException(
            "JointCovWriter requires joint association test results");
    }

    const auto& add = results.additive;
    const auto& dom = *results.dominance;
    const auto& cov = *results.joint_covariance;

    for (size_t i = 0; i < results.freq.size(); ++i)
    {
        const auto row = start + i;
        fmt::format_to(
            std::back_inserter(line_buffer_),
            FMT_COMPILE(
                "{}\t{}\t{}\t{}\t{}\t{:.6g}\t{:.6g}\t{:.6g}\t"
                "{:.6e}\t{:.6e}\t{:.6e}\n"),
            chrom_[row],
            keys_[row],
            pos_[row],
            a1_[row],
            a2_[row],
            results.freq[i],
            add.beta[i],
            dom.beta[i],
            cov.var_a[i],
            cov.var_d[i],
            cov.cov_ad[i]);
    }

    if (line_buffer_.size() >= JOINT_COV_BUFFER_FLUSH_THRESHOLD)
    {
        ofs_.write(
            line_buffer_.data(),
            static_cast<std::streamsize>(line_buffer_.size()));
        line_buffer_.clear();
    }
}

}  // namespace gelex::gwas
