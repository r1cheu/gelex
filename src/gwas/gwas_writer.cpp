#include "gelex/gwas/gwas_writer.h"
#include <cstddef>
#include <ios>
#include <stdexcept>

#include <fmt/compile.h>

namespace gelex::gwas
{

constexpr size_t BUFFER_FLUSH_THRESHOLD = static_cast<const size_t>(64 * 1024);

GwasWriter::GwasWriter(std::string_view out_prefix, AssocMode model)
    : model_(model)
{
    std::string filepath = std::string(out_prefix) + ".gwas.tsv";
    ofs_.open(filepath, std::ios::out | std::ios::binary);

    if (!ofs_.is_open())
    {
        throw std::runtime_error("Failed to open output file: " + filepath);
    }

    line_buffer_.reserve(BUFFER_FLUSH_THRESHOLD);
}

GwasWriter::~GwasWriter()
{
    finalize();
}

auto GwasWriter::write_header() -> void
{
    line_buffer_.clear();

    fmt::format_to(
        std::back_inserter(line_buffer_),
        FMT_COMPILE("CHR\tSNP\tBP\tA1\tA2\tFREQ"));

    switch (model_)
    {
        case AssocMode::A:
            fmt::format_to(
                std::back_inserter(line_buffer_), FMT_COMPILE("\tBETA\tSE"));
            break;
        case AssocMode::D:
            fmt::format_to(
                std::back_inserter(line_buffer_),
                FMT_COMPILE("\tBETA_D\tSE_D"));
            break;
        case AssocMode::AD:
            fmt::format_to(
                std::back_inserter(line_buffer_),
                FMT_COMPILE("\tBETA_A\tSE_A\tBETA_D\tSE_D"));
            break;
    }

    fmt::format_to(std::back_inserter(line_buffer_), FMT_COMPILE("\tSTAT\tP"));
    fmt::format_to(std::back_inserter(line_buffer_), FMT_COMPILE("\tDF\n"));

    ofs_.write(
        line_buffer_.data(), static_cast<std::streamsize>(line_buffer_.size()));
    line_buffer_.clear();
}

auto GwasWriter::write_result(
    const SNPInfo& snp_info,
    const AssociationResult& result) -> void
{
    fmt::format_to(
        std::back_inserter(line_buffer_),
        FMT_COMPILE("{}\t{}\t{}\t{}\t{}\t{:.6g}"),
        snp_info.chrom,
        snp_info.rsid,
        snp_info.bp,
        snp_info.a1,
        snp_info.a2,
        snp_info.freq);

    switch (model_)
    {
        case AssocMode::A:
            fmt::format_to(
                std::back_inserter(line_buffer_),
                FMT_COMPILE("\t{:.6g}\t{:.6g}"),
                result.beta_a,
                result.se_a);
            break;
        case AssocMode::D:
            fmt::format_to(
                std::back_inserter(line_buffer_),
                FMT_COMPILE("\t{:.6g}\t{:.6g}"),
                result.beta_d,
                result.se_d);
            break;
        case AssocMode::AD:
            fmt::format_to(
                std::back_inserter(line_buffer_),
                FMT_COMPILE("\t{:.6g}\t{:.6g}\t{:.6g}\t{:.6g}"),
                result.beta_a,
                result.se_a,
                result.beta_d,
                result.se_d);
            break;
    }

    fmt::format_to(
        std::back_inserter(line_buffer_),
        FMT_COMPILE("\t{:.6g}\t{:.6g}"),
        result.stat,
        result.pvalue);
    fmt::format_to(
        std::back_inserter(line_buffer_), FMT_COMPILE("\t{}\n"), result.df);

    if (line_buffer_.size() >= BUFFER_FLUSH_THRESHOLD)
    {
        ofs_.write(
            line_buffer_.data(),
            static_cast<std::streamsize>(line_buffer_.size()));
        line_buffer_.clear();
    }
}

auto GwasWriter::finalize() -> void
{
    if (ofs_.is_open())
    {
        if (line_buffer_.size() > 0)
        {
            ofs_.write(
                line_buffer_.data(),
                static_cast<std::streamsize>(line_buffer_.size()));
            line_buffer_.clear();
        }
        ofs_.flush();
        ofs_.close();
    }
}

}  // namespace gelex::gwas
