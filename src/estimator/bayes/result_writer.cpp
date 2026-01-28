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
