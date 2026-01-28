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

#include "snp_quant_genetic_writer.h"

#include <format>
#include <fstream>

#include <Eigen/Core>

#include "../src/data/loader/bim_loader.h"
#include "../src/data/parser.h"

namespace gelex
{

using Eigen::Index;
using Eigen::VectorXd;

SnpQuantGeneticWriter::SnpQuantGeneticWriter(
    const MCMCResult& result,
    const std::filesystem::path& bim_file_path)
    : result_(&result), bim_loader_(bim_file_path)
{
}

void SnpQuantGeneticWriter::write(const std::filesystem::path& path) const
{
    if (result_->additive() == nullptr)
    {
        return;
    }

    auto stream = detail::open_file<std::ofstream>(path, std::ios::out);

    write_header(stream);

    // Write original genetic effects for all SNPs
    for (Index i = 0; i < result_->additive()->coeffs.size(); ++i)
    {
        write_snp_row(stream, i);
    }
}

void SnpQuantGeneticWriter::write_header(std::ofstream& stream) const
{
    stream << "Index\tID\tChrom\tPosition\tA1\tA2\tA1Freq\t"
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
    if (snp_index < static_cast<Index>(bim_loader_.size()))
    {
        const auto& snp_info = bim_loader_.info()[snp_index];
        stream << std::format(
            "{}\t{}\t{}\t{}\t{}",
            snp_info.id,
            snp_info.chrom,
            snp_info.pos,
            snp_info.A1,
            snp_info.A2);

        // Calculate A1Freq from genotype mean: mean(X_i) / 2
        double a1_frq = result_->p_freq(snp_index);
        stream << std::format("\t{:.6f}", a1_frq);
        return a1_frq;
    }
    // Placeholder values for missing SNP information
    stream << "\tNA\tNA\tNA\tNA\tNA";  // Chrom, Position, A1, A2, A1Freq
    return std::numeric_limits<double>::quiet_NaN();
}

}  // namespace gelex
