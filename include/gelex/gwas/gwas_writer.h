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

#ifndef GELEX_GWAS_GWAS_WRITER_H
#define GELEX_GWAS_GWAS_WRITER_H

#include <fstream>
#include <string_view>
#include "gelex/types/snp_info.h"

#include <fmt/format.h>

namespace gelex::gwas
{

class GwasWriter
{
   public:
    struct AssocResult
    {
        double freq;
        double beta;
        double se;
        double p_value;
    };
    explicit GwasWriter(std::string_view out_prefix);
    GwasWriter(const GwasWriter&) = delete;
    GwasWriter(GwasWriter&&) = delete;
    GwasWriter& operator=(const GwasWriter&) = delete;
    GwasWriter& operator=(GwasWriter&&) = delete;

    ~GwasWriter();

    auto write_header() -> void;
    auto write_result(const SnpMeta& snp_meta, AssocResult result) -> void;
    auto finalize() -> void;

   private:
    std::ofstream ofs_;

    fmt::memory_buffer line_buffer_;
};

}  // namespace gelex::gwas

#endif  // GELEX_GWAS_GWAS_WRITER_H
