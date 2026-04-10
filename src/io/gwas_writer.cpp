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
#include <stdexcept>
#include <string>

#include <fmt/compile.h>

namespace gelex::gwas
{

constexpr size_t BUFFER_FLUSH_THRESHOLD = static_cast<size_t>(64 * 1024);

GwasWriter::GwasWriter(
    std::string_view out_prefix,
    const df::DataFrame<std::string>& bim)
    : keys_(bim.index().keys()),
      chrom_(bim["chrom"].as<std::string>()),
      pos_(bim["pos"].as<std::int32_t>()),
      a1_(bim["A1"].as<std::string>()),
      a2_(bim["A2"].as<std::string>())
{
    std::string filepath = std::string(out_prefix) + ".gwas.tsv";
    ofs_.open(filepath, std::ios::out | std::ios::binary);

    if (!ofs_.is_open())
    {
        throw std::runtime_error("Failed to open output file: " + filepath);
    }

    line_buffer_.reserve(BUFFER_FLUSH_THRESHOLD);

    fmt::format_to(
        std::back_inserter(line_buffer_),
        FMT_COMPILE("CHR\tSNP\tBP\tA1\tA2\tA1FREQ\tBETA\tSE\tP\tPVE\n"));
    ofs_.write(
        line_buffer_.data(), static_cast<std::streamsize>(line_buffer_.size()));
    line_buffer_.clear();
}

GwasWriter::~GwasWriter()
{
    if (line_buffer_.size() == 0U)
    {
        return;
    }
    ofs_.write(
        line_buffer_.data(), static_cast<std::streamsize>(line_buffer_.size()));
}

auto GwasWriter::write(std::size_t row, AssocResult result) -> void
{
    fmt::format_to(
        std::back_inserter(line_buffer_),
        FMT_COMPILE("{}\t{}\t{}\t{}\t{}\t{:.6g}\t"),
        chrom_[row],
        keys_[row],
        pos_[row],
        a1_[row],
        a2_[row],
        result.freq);

    fmt::format_to(
        std::back_inserter(line_buffer_),
        FMT_COMPILE("{:.6g}\t{:.6g}\t{:.6e}\t{:.6e}\n"),
        result.beta,
        result.se,
        result.p_value,
        result.pve);

    if (line_buffer_.size() >= BUFFER_FLUSH_THRESHOLD)
    {
        ofs_.write(
            line_buffer_.data(),
            static_cast<std::streamsize>(line_buffer_.size()));
        line_buffer_.clear();
    }
}

}  // namespace gelex::gwas
