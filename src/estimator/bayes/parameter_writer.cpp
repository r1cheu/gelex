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
    write_residual_variance(stream);
    write_additive_variance(stream);
    write_dominant_variance(stream);
}

void ParameterWriter::write_fixed_effects(std::ofstream& stream) const
{
    if (result_->fixed() != nullptr)
    {
        write_summary_statistics(
            stream, result_->fixed()->coeffs, result_->fixed()->coeffs.size());
    }
}

void ParameterWriter::write_random_effects(std::ofstream& stream) const
{
    for (const auto& rand : result_->random())
    {
        write_summary_statistics(stream, rand.coeffs, rand.coeffs.size());
        write_summary_statistics(stream, rand.variance, rand.variance.size());
    }
}

void ParameterWriter::write_residual_variance(std::ofstream& stream) const
{
    write_summary_statistics(
        stream, result_->residual(), result_->residual().size());
}

void ParameterWriter::write_additive_variance(std::ofstream& stream) const
{
    if (result_->additive() != nullptr)
    {
        write_summary_statistics(
            stream,
            result_->additive()->variance,
            result_->additive()->variance.size());
    }
}

void ParameterWriter::write_dominant_variance(std::ofstream& stream) const
{
    if (result_->dominant() != nullptr)
    {
        write_summary_statistics(
            stream,
            result_->dominant()->variance,
            result_->dominant()->variance.size());
    }
}

void ParameterWriter::write_summary_statistics(
    std::ofstream& stream,
    const PosteriorSummary& stats,
    Index n_params)
{
    for (Index i = 0; i < n_params; ++i)
    {
        stream << "\t" << stats.mean(i) << "\t" << stats.stddev(i) << "\t"
               << stats.hpdi_low(i) << "\t" << stats.hpdi_high(i) << "\t"
               << stats.ess(i) << "\t" << stats.rhat(i) << "\n";
    }
}

}  // namespace gelex
