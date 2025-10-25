#include "gelex/estimator/bayes/result_writer.h"

#include <format>
#include <fstream>
#include <memory>
#include <optional>
#include <string>

#include <Eigen/Core>
#include "data/loader.h"

namespace gelex
{

using Eigen::Index;
using Eigen::VectorXd;

MCMCResultWriter::MCMCResultWriter(
    const MCMCResult& result,
    const std::filesystem::path& bim_file_path)
    : result_(result)
{
    auto snp_loader = SnpInfoLoader::create(bim_file_path);
    if (snp_loader)
    {
        snp_info_loader_
            = std::make_unique<SnpInfoLoader>(std::move(*snp_loader));
    }
}

void MCMCResultWriter::save(const std::filesystem::path& prefix) const
{
    auto params_path = prefix;
    params_path.replace_extension("params");
    write_parameter_file(params_path);

    auto snp_path = prefix;
    snp_path.replace_extension(".snp.eff");
    write_snp_effects(snp_path);
}

void MCMCResultWriter::write_parameter_file(
    const std::filesystem::path& path) const
{
    auto stream = *detail::open_file<std::ofstream>(path, std::ios_base::out);

    // Write header
    std::string hpdi_low = get_hpdi_low_label();
    std::string hpdi_high = get_hpdi_high_label();
    stream << std::format(
        "term\tmean\tstddev\t{}%\t{}%\tess\trhat\n", hpdi_low, hpdi_high);

    // Write fixed effects
    if (result_.fixed() != nullptr)
    {
        write_summary_statistics(
            stream, result_.fixed()->coeffs, result_.fixed()->coeffs.size());
    }

    // Write random effects
    for (const auto& rand : result_.random())
    {
        write_summary_statistics(stream, rand.coeffs, rand.coeffs.size());
        write_summary_statistics(stream, rand.variance, rand.variance.size());
    }

    // Write residual variance
    write_summary_statistics(
        stream, result_.residual(), result_.residual().size());

    // Write additive variance
    if (result_.additive() != nullptr)
    {
        write_summary_statistics(
            stream,
            result_.additive()->variance,
            result_.additive()->variance.size());
    }

    // Write dominant variance
    if (result_.dominant() != nullptr)
    {
        write_summary_statistics(
            stream,
            result_.dominant()->variance,
            result_.dominant()->variance.size());
    }
}

void MCMCResultWriter::write_summary_statistics(
    std::ofstream& stream,
    const PosteriorSummary& stats,
    Index n_params) const
{
    for (Index i = 0; i < n_params; ++i)
    {
        stream << "\t" << stats.mean(i) << "\t" << stats.stddev(i) << "\t"
               << stats.hpdi_low(i) << "\t" << stats.hpdi_high(i) << "\t"
               << stats.ess(i) << "\t" << stats.rhat(i) << "\n";
    }
}

void MCMCResultWriter::write_snp_effects(
    const std::filesystem::path& path) const
{
    if (result_.additive() == nullptr)
    {
        return;
    }

    auto stream = *detail::open_file<std::ofstream>(path, std::ios_base::out);

    // Determine if we have per-component probabilities to write
    const Eigen::Index n_components = result_.snp_tracker()->comp_probs.cols();

    // Write dynamic header based on whether dominance effects exist and
    // component probs
    stream << "Index\tID\tChrom\tPosition\tA1\tA2\tA1Frq\tAdd\tAddSE\tAddPVE";
    if (n_components > 2)
    {
        for (Eigen::Index comp = 0; comp < n_components; ++comp)
        {
            stream << "\tVg" << comp << "";
        }
        stream << "\tPIP";
    }
    else
    {
        stream << "\tPIP";
    }

    if (result_.dominant() != nullptr)
    {
        stream << "tDomEff\tDomSE\tDomPVE\td / a";
    }

    stream << "\n";

    // Write SNP effects for all SNPs
    for (Index i = 0; i < result_.additive()->coeffs.size(); ++i)
    {
        stream << std::format("{}\t", i + 1);  // Index

        // SNP name and basic information
        if (snp_info_loader_
            && i < static_cast<Index>(snp_info_loader_->size()))
        {
            const auto& snp_info = (*snp_info_loader_)[i];
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
            else if (result_.additive_means_.size() > i)
            {
                // Calculate A1Frq from genotype mean: mean(X_i) / 2
                double a1_frq = result_.additive_means_(i) / 2.0;
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

        // Additive effect statistics
        stream << std::format(
            "\t{:.6f}\t{:.6f}",
            result_.additive()->coeffs.mean(i),
            result_.additive()->coeffs.stddev(i));

        // Additive PVE statistics
        if (result_.additive()->pve.size() > i)
        {
            stream << std::format("\t{:.6e}", result_.additive()->pve.mean(i));
        }
        else
        {
            stream << "\t0.0";  // Placeholder for PVE
        }

        // Per-component posterior probabilities (if available)
        if (i < result_.snp_tracker()->comp_probs.rows())
        {
            for (Eigen::Index comp = 0; comp < n_components; ++comp)
            {
                stream << std::format(
                    "\t{:.6f}", (result_.snp_tracker()->comp_probs)(i, comp));
            }
        }

        // Posterior inclusion probability (if available)
        if ((result_.snp_tracker() != nullptr)
            && result_.snp_tracker()->pip.size() > i)
        {
            stream << std::format("\t{:.6f}", result_.snp_tracker()->pip(i));
        }
        else
        {
            stream << "\t1.0";  // Default PIP when not tracking
        }

        // Dominance effect statistics (if they exist)
        if ((result_.dominant() != nullptr)
            && i < result_.dominant()->coeffs.size())
        {
            stream << std::format(
                "\t{:.6f}\t{:.6f}",
                result_.dominant()->coeffs.mean(i),
                result_.dominant()->coeffs.stddev(i));

            // Dominant PVE statistics
            if (result_.dominant()->pve.size() > i)
            {
                stream << std::format(
                    "\t{:.6e}", result_.dominant()->pve.mean(i));
            }
            else
            {
                stream << "\t0.0";  // Placeholder for PVE
            }

            // Ratio (d/a) statistics
            if (result_.dominant()->ratios.size() > i)
            {
                stream << std::format(
                    "\t{:.6f}", result_.dominant()->ratios.mean(i));
            }
            else
            {
                stream << "\t0.0";  // Placeholder for ratio
            }
        }

        stream << "\n";
    }
}

std::string MCMCResultWriter::get_hpdi_low_label() const
{
    return std::to_string(
        static_cast<int>(
            std::round(100 * (1 - result_.residual().mean(0)) / 2)));
}

std::string MCMCResultWriter::get_hpdi_high_label() const
{
    return std::to_string(
        static_cast<int>(
            std::round(100 * (1 + result_.residual().mean(0)) / 2)));
}

}  // namespace gelex
