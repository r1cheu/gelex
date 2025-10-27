#include "snp_effects_writer.h"

#include <format>
#include <fstream>
#include <memory>

#include <Eigen/Core>
#include "data/loader.h"

namespace gelex
{

using Eigen::Index;
using Eigen::VectorXd;

SnpEffectsWriter::SnpEffectsWriter(
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

void SnpEffectsWriter::write(const std::filesystem::path& path) const
{
    if (result_->additive() == nullptr)
    {
        return;
    }

    auto stream = *detail::open_file<std::ofstream>(path, std::ios_base::out);

    write_header(stream);

    // Write SNP effects for all SNPs
    for (Index i = 0; i < result_->additive()->coeffs.size(); ++i)
    {
        write_snp_row(stream, i);
    }
}

void SnpEffectsWriter::write_header(std::ofstream& stream) const
{
    // Determine if we have per-component probabilities to write
    const Index n_components = result_->snp_tracker()->comp_probs.cols();

    // Write dynamic header based on whether dominance effects exist and
    // component probs
    stream << "Index\tID\tChrom\tPosition\tA1\tA2\tA1Frq\tAdd\tAddSE\tAddPVE";
    if (n_components > 2)
    {
        for (Index comp = 0; comp < n_components; ++comp)
        {
            stream << "\tpi_" << comp << "";
        }
        stream << "\tPIP";
    }
    else
    {
        stream << "\tPIP";
    }

    if (has_dominant_effects())
    {
        stream << "tDomEff\tDomSE\tDomPVE\td / a";
    }

    stream << "\n";
}

void SnpEffectsWriter::write_snp_row(std::ofstream& stream, Index snp_index)
    const
{
    stream << std::format("{}\t", snp_index + 1);  // Index

    write_snp_basic_info(stream, snp_index);
    write_additive_effects(stream, snp_index);
    write_component_probabilities(stream, snp_index);
    write_pip(stream, snp_index);
    write_dominant_effects(stream, snp_index);

    stream << "\n";
}

void SnpEffectsWriter::write_snp_basic_info(
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

        // Use allele frequency if available, otherwise calculate from
        // genotype means
        if (snp_info.a1_frq)
        {
            stream << std::format("\t{}", *snp_info.a1_frq);
        }
        else if (result_->additive_means_.size() > snp_index)
        {
            // Calculate A1Frq from genotype mean: mean(X_i) / 2
            double a1_frq = result_->additive_means_(snp_index) / 2.0;
            stream << std::format("\t{:.6f}", a1_frq);
        }
        else
        {
            stream << "\t0.5";  // Placeholder for A1Frq
        }
    }
    else
    {
        // Placeholder values for missing SNP information
        stream << "\t1\t1000\tA\tC\t0.5";  // Chrom, Position, A1, A2, A1Frq
    }
}

void SnpEffectsWriter::write_additive_effects(
    std::ofstream& stream,
    Index snp_index) const
{
    // Additive effect statistics
    stream << std::format(
        "\t{:.6f}\t{:.6f}",
        result_->additive()->coeffs.mean(snp_index),
        result_->additive()->coeffs.stddev(snp_index));

    // Additive PVE statistics
    if (result_->additive()->pve.size() > snp_index)
    {
        stream << std::format(
            "\t{:.6e}", result_->additive()->pve.mean(snp_index));
    }
    else
    {
        stream << "\t0.0";  // Placeholder for PVE
    }
}

void SnpEffectsWriter::write_component_probabilities(
    std::ofstream& stream,
    Index snp_index) const
{
    // Per-component posterior probabilities (if available)
    if (has_component_probabilities()
        && snp_index < result_->snp_tracker()->comp_probs.rows())
    {
        const Index n_components = result_->snp_tracker()->comp_probs.cols();
        for (Index comp = 0; comp < n_components; ++comp)
        {
            stream << std::format(
                "\t{:.6f}",
                (result_->snp_tracker()->comp_probs)(snp_index, comp));
        }
    }
}

void SnpEffectsWriter::write_pip(std::ofstream& stream, Index snp_index) const
{
    // Posterior inclusion probability (if available)
    if ((result_->snp_tracker() != nullptr)
        && result_->snp_tracker()->pip.size() > snp_index)
    {
        stream << std::format(
            "\t{:.6f}", result_->snp_tracker()->pip(snp_index));
    }
    else
    {
        stream << "\t1.0";  // Default PIP when not tracking
    }
}

void SnpEffectsWriter::write_dominant_effects(
    std::ofstream& stream,
    Index snp_index) const
{
    // Dominance effect statistics (if they exist)
    if (has_dominant_effects()
        && snp_index < result_->dominant()->coeffs.size())
    {
        stream << std::format(
            "\t{:.6f}\t{:.6f}",
            result_->dominant()->coeffs.mean(snp_index),
            result_->dominant()->coeffs.stddev(snp_index));

        // Dominant PVE statistics
        if (result_->dominant()->pve.size() > snp_index)
        {
            stream << std::format(
                "\t{:.6e}", result_->dominant()->pve.mean(snp_index));
        }
        else
        {
            stream << "\t0.0";  // Placeholder for PVE
        }

        // Ratio (d/a) statistics
        if (result_->dominant()->ratios.size() > snp_index)
        {
            stream << std::format(
                "\t{:.6f}", result_->dominant()->ratios.mean(snp_index));
        }
        else
        {
            stream << "\t0.0";  // Placeholder for ratio
        }
    }
}

double SnpEffectsWriter::get_allele_frequency(Index snp_index) const
{
    if (snp_info_loader_
        && snp_index < static_cast<Index>(snp_info_loader_->size()))
    {
        const auto& snp_info = (*snp_info_loader_)[snp_index];
        if (snp_info.a1_frq)
        {
            return *snp_info.a1_frq;
        }
    }

    if (result_->additive_means_.size() > snp_index)
    {
        return result_->additive_means_(snp_index) / 2.0;
    }

    return 0.5;  // Default allele frequency
}

bool SnpEffectsWriter::has_component_probabilities() const
{
    return (result_->snp_tracker() != nullptr)
           && result_->snp_tracker()->comp_probs.rows() > 0;
}

bool SnpEffectsWriter::has_dominant_effects() const
{
    return result_->dominant() != nullptr;
}

}  // namespace gelex
