#ifndef GELEX_GWAS_GWAS_WRITER_H
#define GELEX_GWAS_GWAS_WRITER_H

#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>

#include "gelex/gwas/association_test.h"
#include "gelex/gwas/snp_encoder.h"

namespace gelex::gwas
{

struct SNPInfo
{
    std::string chrom;
    std::string rsid;
    int64_t bp{};
    std::string a1;  // effect allele
    std::string a2;  // other allele
    double freq{};   // frequency of A1
    int n{};         // sample size
};

class GwasWriter
{
   public:
    GwasWriter(
        std::string_view out_prefix,
        GwasModel model,
        TestType test_type);

    auto write_header() -> void;
    auto write_result(const SNPInfo& snp_info, const AssociationResult& result)
        -> void;
    auto finalize() -> void;

   private:
    std::ofstream ofs_;
    GwasModel model_;
    TestType test_type_;
};

}  // namespace gelex::gwas

#endif  // GELEX_GWAS_GWAS_WRITER_H
