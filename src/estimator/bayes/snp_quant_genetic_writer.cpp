#include "snp_quant_genetic_writer.h"

#include <format>
#include <fstream>
#include <memory>

#include <Eigen/Core>
#include "data/loader.h"

namespace gelex
{

using Eigen::Index;
using Eigen::VectorXd;

SnpQuantGeneticWriter::SnpQuantGeneticWriter(
    const MCMCResult& result,
    const std::filesystem::path& bim_file_path)
    : result_(&result)
{
    auto snp_loader = SnpInfoLoader::create(bim_file_path);
    if (snp_loader)
    {
        snp_info_loader_
            = std::make_unique<SnpInfoLoader>(std::move(*snp_loader));
    }
}

void SnpQuantGeneticWriter::write(const std::filesystem::path& path) const
{
    if (result_->additive() == nullptr)
    {
        return;
    }

    auto stream = *detail::open_file<std::ofstream>(path, std::ios_base::out);

    write_header(stream);

    // Write original genetic effects for all SNPs
    for (Index i = 0; i < result_->additive()->coeffs.size(); ++i)
    {
        write_snp_row(stream, i);
    }
}

void SnpQuantGeneticWriter::write_header(std::ofstream& stream) const
{
    stream << "Index\tID\tChrom\tPosition\tA1\tA2\tA1Frq\t"
           << "a\t"
           << "d\t"
           << "d/|a|\n";
}

void SnpQuantGeneticWriter::write_snp_row(
    std::ofstream& stream,
    Index snp_index) const
{
    stream << std::format("{}\t", snp_index + 1);  // Index

    const double p_freq = write_snp_basic_info(stream, snp_index);
    const double q_freq = 1.0 - p_freq;
    const double delta = result_->dominant() != nullptr
                             ? result_->dominant()->coeffs.mean(snp_index)
                             : 0.0;
    const double alpha = result_->additive()->coeffs.mean(snp_index);

    const double d_stddev = 2 * p_freq * q_freq;
    const double a_stddev = std::sqrt(d_stddev);

    const double d = delta / d_stddev;
    const double a = (alpha / a_stddev) + ((p_freq - q_freq) * d);

    const double d_over_a = (a != 0.0) ? (d / std::abs(a)) : 0.0;

    stream << std::format("\t{:.6f}\t{:.6f}\t{:.6f}", a, d, d_over_a);

    stream << "\n";
}

double SnpQuantGeneticWriter::write_snp_basic_info(
    std::ofstream& stream,
    Index snp_index) const
{
    // SNP name and basic information
    if (snp_info_loader_
        && snp_index < static_cast<Index>(snp_info_loader_->size()))
    {
        const auto& snp_info = (*snp_info_loader_)[snp_index];
        stream << std::format(
            "{}\t{}\t{}\t{}\t{}",
            snp_info.id,
            snp_info.chrom,
            snp_info.position,
            snp_info.a1,
            snp_info.a2);

        // Calculate A1Frq from genotype mean: mean(X_i) / 2
        double a1_frq = result_->p_freq(snp_index);
        stream << std::format("\t{:.6f}", a1_frq);
        return a1_frq;
    }
    // Placeholder values for missing SNP information
    stream << "\tNA\tNA\tNA\tNA\tNA";  // Chrom, Position, A1, A2, A1Frq
    return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace gelex
