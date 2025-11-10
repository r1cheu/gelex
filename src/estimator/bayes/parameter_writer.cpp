#include "parameter_writer.h"

#include <fstream>

#include <Eigen/Core>

#include "data/loader.h"

namespace gelex
{

using Eigen::Index;
using Eigen::VectorXd;

ParameterWriter::ParameterWriter(const MCMCResult& result) : result_(&result) {}

void ParameterWriter::write(const std::filesystem::path& path) const
{
    auto stream = *detail::open_file<std::ofstream>(path, std::ios_base::out);

    // Write header
    stream << "term\tmean\tstddev\t5%\t95%\tess\trhat\n";

    // Write different parameter types
    write_fixed_effects(stream);
    write_random_effects(stream);
    write_additive_effect(stream);
    write_dominant_effect(stream);
    write_residual_variance(stream);
}

void ParameterWriter::write_fixed_effects(std::ofstream& stream) const
{
    if (result_->fixed() != nullptr)
    {
        std::vector<std::string> terms;
        terms.reserve(result_->fixed()->coeffs.size());
        for (int i = 0; i < result_->fixed()->coeffs.size(); ++i)
        {
            terms.emplace_back("Intercept");
        }

        write_summary_statistics(
            terms,
            stream,
            result_->fixed()->coeffs,
            result_->fixed()->coeffs.size());
    }
}

void ParameterWriter::write_random_effects(std::ofstream& stream) const
{
    for (const auto& rand : result_->random())
    {
        write_summary_statistics(
            std::vector<std::string>(rand.coeffs.size()),
            stream,
            rand.coeffs,
            rand.coeffs.size());
        write_summary_statistics(
            std::vector<std::string>(rand.variance.size()),
            stream,
            rand.variance,
            rand.variance.size());
    }
}

void ParameterWriter::write_residual_variance(std::ofstream& stream) const
{
    write_summary_statistics(
        std::vector<std::string>{"σ²_e"},
        stream,
        result_->residual(),
        result_->residual().size());
}

void ParameterWriter::write_additive_effect(std::ofstream& stream) const
{
    if (result_->additive() != nullptr)
    {
        write_summary_statistics(
            std::vector<std::string>{"σ²_add"},
            stream,
            result_->additive()->variance,
            result_->additive()->variance.size());
        const auto& mixture_proportion
            = result_->additive()->mixture_proportion;
        std::vector<std::string> proportion_terms;
        proportion_terms.reserve(mixture_proportion.size());
        for (Index i = 0; i < mixture_proportion.size(); ++i)
        {
            proportion_terms.emplace_back("π[" + std::to_string(1) + "]");
        }
        write_summary_statistics(
            std::vector<std::string>{"h²"},
            stream,
            result_->additive()->heritability,
            result_->additive()->heritability.size());
        write_summary_statistics(
            proportion_terms,
            stream,
            mixture_proportion,
            mixture_proportion.size());
    }
}

void ParameterWriter::write_dominant_effect(std::ofstream& stream) const
{
    if (result_->dominant() != nullptr)
    {
        write_summary_statistics(
            std::vector<std::string>{"σ²_dom"},
            stream,
            result_->dominant()->variance,
            result_->dominant()->variance.size());
        const auto& mixture_proportion
            = result_->dominant()->mixture_proportion;
        std::vector<std::string> proportion_terms;
        proportion_terms.reserve(mixture_proportion.size());
        for (Index i = 0; i < mixture_proportion.size(); ++i)
        {
            proportion_terms.emplace_back("π[" + std::to_string(1) + "]");
        }
        write_summary_statistics(
            std::vector<std::string>{"δ²"},
            stream,
            result_->dominant()->heritability,
            result_->dominant()->heritability.size());
        write_summary_statistics(
            proportion_terms,
            stream,
            mixture_proportion,
            mixture_proportion.size());
    }
}

void ParameterWriter::write_summary_statistics(
    std::span<const std::string> terms,
    std::ofstream& stream,
    const PosteriorSummary& stats,
    Index n_params)
{
    for (Index i = 0; i < n_params; ++i)
    {
        stream << terms[i] << "\t" << stats.mean(i) << "\t" << stats.stddev(i)
               << "\t" << stats.hpdi_low(i) << "\t" << stats.hpdi_high(i)
               << "\t" << stats.ess(i) << "\t" << stats.rhat(i) << "\n";
    }
}

}  // namespace gelex
