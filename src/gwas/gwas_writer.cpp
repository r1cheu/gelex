#include "gelex/gwas/gwas_writer.h"

#include <fmt/format.h>
#include <iomanip>
#include <stdexcept>

namespace gelex::gwas
{

GwasWriter::GwasWriter(
    std::string_view out_prefix,
    GwasModel model,
    TestType test_type)
    : model_(model), test_type_(test_type)
{
    std::string filepath = std::string(out_prefix) + ".gwas.tsv";
    ofs_.open(filepath);
    if (!ofs_.is_open())
    {
        throw std::runtime_error("Failed to open output file: " + filepath);
    }
}

auto GwasWriter::write_header() -> void
{
    ofs_ << "CHR\tSNP\tBP\tA1\tA2\tFREQ";

    switch (model_)
    {
        case GwasModel::Additive:
            ofs_ << "\tBETA\tSE";
            break;
        case GwasModel::Dominance:
            ofs_ << "\tBETA_D\tSE_D";
            break;
        case GwasModel::AdditiveDominance:
            ofs_ << "\tBETA_A\tSE_A\tBETA_D\tSE_D";
            break;
    }

    ofs_ << "\tSTAT\tP";

    // Add separate p-values if test type is separate
    if (model_ == GwasModel::AdditiveDominance
        && test_type_ == TestType::Separate)
    {
        ofs_ << "\tP_A\tP_D";
    }

    ofs_ << "\tDF\tN\n";
}

auto GwasWriter::write_result(
    const SNPInfo& snp_info,
    const AssociationResult& result) -> void
{
    // CHR, SNP, BP, A1, A2, FREQ
    ofs_ << snp_info.chrom << "\t" << snp_info.rsid << "\t" << snp_info.bp
         << "\t" << snp_info.a1 << "\t" << snp_info.a2 << "\t";

    ofs_ << fmt::format("{:.6g}", snp_info.freq);

    // Effect estimates based on model
    switch (model_)
    {
        case GwasModel::Additive:
            ofs_ << "\t" << fmt::format("{:.6g}", result.beta_a) << "\t"
                 << fmt::format("{:.6g}", result.se_a);
            break;
        case GwasModel::Dominance:
            ofs_ << "\t" << fmt::format("{:.6g}", result.beta_d) << "\t"
                 << fmt::format("{:.6g}", result.se_d);
            break;
        case GwasModel::AdditiveDominance:
            ofs_ << "\t" << fmt::format("{:.6g}", result.beta_a) << "\t"
                 << fmt::format("{:.6g}", result.se_a) << "\t"
                 << fmt::format("{:.6g}", result.beta_d) << "\t"
                 << fmt::format("{:.6g}", result.se_d);
            break;
    }

    // STAT, P
    ofs_ << "\t" << fmt::format("{:.6g}", result.stat) << "\t"
         << fmt::format("{:.6g}", result.pvalue);

    // Separate p-values if applicable
    if (model_ == GwasModel::AdditiveDominance
        && test_type_ == TestType::Separate)
    {
        ofs_ << "\t" << fmt::format("{:.6g}", result.pvalue_a) << "\t"
             << fmt::format("{:.6g}", result.pvalue_d);
    }

    // DF, N
    ofs_ << "\t" << result.df << "\t" << snp_info.n << "\n";
}

auto GwasWriter::finalize() -> void
{
    ofs_.flush();
    ofs_.close();
}

}  // namespace gelex::gwas
