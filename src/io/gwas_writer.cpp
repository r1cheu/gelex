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

#include "gelex/io/gwas_writer.h"

#include <cstddef>
#include <cstdint>
#include <ios>
#include <string>

#include <fmt/compile.h>

#include "gelex/algo/gwas/assoc_tester.h"

namespace gelex::gwas
{

constexpr size_t BUFFER_FLUSH_THRESHOLD = static_cast<size_t>(64 * 1024);

GwasWriter::GwasWriter(
    std::string_view out_prefix,
    const dataframe::DataFrame<std::string>& bim,
    AssocTestType test_type)
    : test_type_(test_type),
      keys_(bim.index().keys()),
      chrom_(bim["chrom"].as<std::string>()),
      pos_(bim["pos"].as<std::int32_t>()),
      a1_(bim["A1"].as<std::string>()),
      a2_(bim["A2"].as<std::string>()),
      ofs_(
          std::string(out_prefix) + ".gwas.tsv",
          std::ios::out | std::ios::binary)
{
    line_buffer_.reserve(BUFFER_FLUSH_THRESHOLD);

    switch (test_type_)
    {
        case AssocTestType::Single:
            fmt::format_to(
                std::back_inserter(line_buffer_),
                FMT_COMPILE(
                    "CHR\tSNP\tBP\tA1\tA2\tA1FREQ\tBETA\tSE\tP\tPVE\n"));
            break;
        case AssocTestType::Joint:
            fmt::format_to(
                std::back_inserter(line_buffer_),
                FMT_COMPILE(
                    "CHR\tSNP\tBP\tA1\tA2\tA1FREQ\t"
                    "BETA_A\tSE_A\tP_A\tPVE_A\t"
                    "BETA_D\tSE_D\tP_D\tPVE_D\tPVE\n"));
            break;
    }

    ofs_.write(
        line_buffer_.data(), static_cast<std::streamsize>(line_buffer_.size()));
    line_buffer_.clear();
}

GwasWriter::~GwasWriter() noexcept
{
    if (line_buffer_.size() == 0U)
    {
        return;
    }
    try
    {
        ofs_.write(
            line_buffer_.data(),
            static_cast<std::streamsize>(line_buffer_.size()));
    }
    catch (...)  // NOLINT(bugprone-empty-catch): dtor must be noexcept
    {
        ofs_.setstate(std::ios::failbit);
    }
}

auto GwasWriter::write(std::size_t start, const TestResults& results) -> void
{
    const auto& add = results.additive;
    const bool is_joint = test_type_ == AssocTestType::Joint;

    for (size_t i = 0; i < results.freq.size(); ++i)
    {
        const auto row = start + i;
        fmt::format_to(
            std::back_inserter(line_buffer_),
            FMT_COMPILE("{}\t{}\t{}\t{}\t{}\t{:.6g}\t"),
            chrom_[row],
            keys_[row],
            pos_[row],
            a1_[row],
            a2_[row],
            results.freq[i]);

        fmt::format_to(
            std::back_inserter(line_buffer_),
            FMT_COMPILE("{:.6g}\t{:.6g}\t{:.6e}\t{:.6e}"),
            add.beta[i],
            add.se[i],
            add.p[i],
            add.pve[i]);

        if (is_joint)
        {
            const auto& dom = *results.dominance;
            const auto& total_pve = *results.total_pve;
            fmt::format_to(
                std::back_inserter(line_buffer_),
                FMT_COMPILE("\t{:.6g}\t{:.6g}\t{:.6e}\t{:.6e}\t{:.6e}"),
                dom.beta[i],
                dom.se[i],
                dom.p[i],
                dom.pve[i],
                total_pve[i]);
        }

        line_buffer_.push_back('\n');
    }

    if (line_buffer_.size() >= BUFFER_FLUSH_THRESHOLD)
    {
        ofs_.write(
            line_buffer_.data(),
            static_cast<std::streamsize>(line_buffer_.size()));
        line_buffer_.clear();
    }
}

}  // namespace gelex::gwas
