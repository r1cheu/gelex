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

#include "gelex/pipeline/report/result_writer.h"

#include <filesystem>

#include "gelex/pipeline/report/parameter_writer.h"
#include "gelex/pipeline/report/snp_effects_writer.h"
#include "gelex/types/mcmc_results.h"

namespace gelex
{

MCMCResultWriter::MCMCResultWriter(
    const MCMCResult& result,
    const std::filesystem::path& bim_file_path)
    : result_(&result), bim_file_path_(bim_file_path)
{
}

auto MCMCResultWriter::save(const std::filesystem::path& prefix) const -> void
{
    auto params_path = prefix;
    params_path.replace_extension("params");
    ParameterWriter parameter_writer(*result_, params_path);
    parameter_writer.write();

    if (result_->additive() != nullptr)
    {
        auto snp_path = prefix;
        snp_path.replace_extension(".snp.eff");
        SnpEffectsWriter snp_effects_writer(*result_, bim_file_path_, snp_path);
        snp_effects_writer.write();
    }
}

}  // namespace gelex
