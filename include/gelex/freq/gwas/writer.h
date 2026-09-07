// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_GWAS_WRITER_H_
#define GELEX_FREQ_GWAS_WRITER_H_

#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <span>
#include <string>
#include <string_view>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/freq/gwas/assoc_type.h"
#include "gelex/io/detail/atomic_output_stream.h"

namespace gelex
{
struct TestResults;
}

namespace gelex
{

class GwasWriter
{
   public:
    GwasWriter(
        std::string_view out_prefix,
        const DataFrame<std::string>& bim,
        AssocType test_type = AssocType::Single);
    GwasWriter(const GwasWriter&) = delete;
    GwasWriter(GwasWriter&&) = delete;
    GwasWriter& operator=(const GwasWriter&) = delete;
    GwasWriter& operator=(GwasWriter&&) = delete;

    ~GwasWriter() noexcept;

    auto write(std::size_t start, const TestResults& results) -> void;

   private:
    AssocType test_type_;
    fmt::memory_buffer line_buffer_;

    std::span<const std::string> keys_;
    std::span<const std::string> chrom_;
    std::span<const std::int32_t> pos_;
    std::span<const std::string> a1_;
    std::span<const std::string> a2_;

    detail::AtomicOutputStream ofs_;
};

}  // namespace gelex

#endif  // GELEX_FREQ_GWAS_WRITER_H_
