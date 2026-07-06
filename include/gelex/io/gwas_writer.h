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

#ifndef GELEX_IO_GWAS_WRITER_H_
#define GELEX_IO_GWAS_WRITER_H_

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>

#include <fmt/format.h>

#include "gelex/algo/gwas/assoc_type.h"
#include "gelex/data/dataframe/dataframe.h"
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

#endif  // GELEX_IO_GWAS_WRITER_H_
