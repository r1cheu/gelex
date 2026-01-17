#ifndef GELEX_GWAS_GWAS_WRITER_H
#define GELEX_GWAS_GWAS_WRITER_H

#include <fstream>
#include <string>
#include <string_view>

#include <fmt/format.h>

#include "gelex/gwas/association_test.h"
#include "gelex/gwas/snp_encoder.h"

namespace gelex::gwas
{

struct SNPInfo
{
    std::string chrom;
    std::string rsid;
    int64_t bp{};
    std::string a1;
    std::string a2;
    double freq{};
    int n{};
};

class GwasWriter
{
   public:
    GwasWriter(
        std::string_view out_prefix,
        GwasModel model,
        TestType test_type);
    GwasWriter(const GwasWriter&) = delete;
    GwasWriter(GwasWriter&&) = delete;
    GwasWriter& operator=(const GwasWriter&) = delete;
    GwasWriter& operator=(GwasWriter&&) = delete;

    ~GwasWriter();

    auto write_header() -> void;
    auto write_result(const SNPInfo& snp_info, const AssociationResult& result)
        -> void;
    auto finalize() -> void;

   private:
    std::ofstream ofs_;
    GwasModel model_;
    TestType test_type_;

    fmt::memory_buffer line_buffer_;
};

}  // namespace gelex::gwas

#endif  // GELEX_GWAS_GWAS_WRITER_H
