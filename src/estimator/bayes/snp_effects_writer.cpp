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
    Index n_alpha_components = 0;
    Index n_dominant_components = 0;
    if (const auto* additive = result_->additive();
        (additive != nullptr) && additive->comp_probs.cols() > 0)
    {
        n_alpha_components = additive->comp_probs.cols();
    }
    if (const auto* dominant = result_->dominant();
        (dominant != nullptr) && dominant->comp_probs.cols() > 0)
    {
        n_dominant_components = dominant->comp_probs.cols();
    }

    // Write dynamic header based on whether dominance effects exist and
    // component probs
    stream << "Index\tID\tChrom\tPosition\tA1\tA2\tA1Frq\tAdd\tAddSE\tAddPVE";

    // Write additive component probability columns for any number of components
    if (n_alpha_components > 2)
    {
        for (Index comp = 0; comp < n_alpha_components; ++comp)
        {
            stream << "\tpi_" << comp << "";
        }
    }
    stream << "\tPIP";

    if (result_->dominant() != nullptr)
    {
        stream << "\tDomEff\tDomSE\tDomPVE";
        // Write dominant component probability columns for any number of
        // components
        if (n_dominant_components > 2)
        {
            for (Index comp = 0; comp < n_dominant_components; ++comp)
            {
                stream << "\tpi_" << comp << "";
            }
        }
        stream << "\tPIP";
    }

    stream << "\n";
}

void SnpEffectsWriter::write_snp_row(std::ofstream& stream, Index snp_index)
    const
{
    stream << std::format("{}\t", snp_index + 1);  // Index

    write_snp_basic_info(stream, snp_index);
    write_additive_effects(stream, snp_index);
    write_add_component_probabilities(stream, snp_index);
    write_add_pip(stream, snp_index);
    write_dominant_effects(stream, snp_index);
    write_dom_component_probabilities(stream, snp_index);
    write_dom_pip(stream, snp_index);

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

        if (result_->p_freq.size() > snp_index)
        {
            // Calculate A1Frq from genotype mean: mean(X_i) / 2
            double a1_frq = result_->p_freq(snp_index);
            stream << std::format("\t{:.6f}", a1_frq);
        }
        else
        {
            stream << "\tNA";  // Placeholder for A1Frq
        }
    }
    else
    {
        // Placeholder values for missing SNP information
        stream << "\tNA\tNA\tNA\tNA\tNA";  // Chrom, Position, A1, A2, A1Frq
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

void SnpEffectsWriter::write_add_component_probabilities(
    std::ofstream& stream,
    Index snp_index) const
{
    // Per-component posterior probabilities (if available)
    if (const auto* additive = result_->additive();
        (additive != nullptr) && additive->comp_probs.cols() > 2
        && snp_index < additive->comp_probs.rows())
    {
        const Index n_components = additive->comp_probs.cols();

        for (Index comp = 0; comp < n_components; ++comp)
        {
            stream << std::format(
                "\t{:.6f}", additive->comp_probs(snp_index, comp));
        }
    }
}

void SnpEffectsWriter::write_add_pip(std::ofstream& stream, Index snp_index)
    const
{
    // Posterior inclusion probability (if available)
    if (const auto* additive = result_->additive();
        (additive != nullptr) && additive->pip.size() > snp_index)
    {
        stream << std::format("\t{:.6f}", additive->pip(snp_index));
    }
    else
    {
        stream << "\t1.0";  // Default PIP when not tracking
    }
}

void SnpEffectsWriter::write_dom_component_probabilities(
    std::ofstream& stream,
    Index snp_index) const
{
    // Per-component posterior probabilities (if available)

    if (const auto* dominant = result_->dominant();
        (dominant != nullptr) && dominant->comp_probs.cols() > 2
        && snp_index < dominant->comp_probs.rows())
    {
        const Index n_components = dominant->comp_probs.cols();
        for (Index comp = 0; comp < n_components; ++comp)
        {
            stream << std::format(
                "\t{:.6f}", dominant->comp_probs(snp_index, comp));
        }
    }
}

void SnpEffectsWriter::write_dom_pip(std::ofstream& stream, Index snp_index)
    const
{
    if (const auto* dominant = result_->dominant();
        (dominant != nullptr) && dominant->pip.size() > snp_index)
    {
        stream << std::format("\t{:.6f}", dominant->pip(snp_index));
    }
    else if (const auto* dominant = result_->dominant(); dominant != nullptr)
    {
        stream << "\t1.0";  // Default PIP when not tracking
    }
}

void SnpEffectsWriter::write_dominant_effects(
    std::ofstream& stream,
    Index snp_index) const
{
    // Dominance effect statistics (if they exist)
    if (const auto* dominant = result_->dominant();
        (dominant != nullptr) && snp_index < result_->dominant()->coeffs.size())
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
    }
}
}  // namespace gelex
