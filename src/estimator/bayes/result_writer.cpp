#include "gelex/estimator/bayes/result_writer.h"

#include <filesystem>

namespace gelex
{

MCMCResultWriter::MCMCResultWriter(
    const MCMCResult& result,
    const std::filesystem::path& bim_file_path)
    : parameter_writer_(result),
      snp_effects_writer_(result, bim_file_path),
      snp_quant_genetic_writer_(result, bim_file_path)
{
}

void MCMCResultWriter::save(const std::filesystem::path& prefix) const
{
    auto params_path = prefix;
    params_path.replace_extension("params");
    parameter_writer_.write(params_path);

    auto snp_path = prefix;
    snp_path.replace_extension(".snp.eff");
    snp_effects_writer_.write(snp_path);

    auto quant_path = prefix;
    quant_path.replace_extension(".snp.quant.eff");
    snp_quant_genetic_writer_.write(quant_path);
}

}  // namespace gelex
