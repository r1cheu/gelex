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
