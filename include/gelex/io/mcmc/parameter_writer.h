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

#ifndef GELEX_IO_MCMC_PARAMETER_WRITER_H_
#define GELEX_IO_MCMC_PARAMETER_WRITER_H_

#include <filesystem>
#include <memory>
#include <span>
#include <string>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/result.h"

namespace gelex::io::detail
{
class TextWriter;
}

namespace gelex
{

class ParameterWriter
{
   public:
    ParameterWriter(
        const mcmc::Result& result,
        const std::filesystem::path& output_path);
    ~ParameterWriter();
    ParameterWriter(const ParameterWriter&) = delete;
    ParameterWriter& operator=(const ParameterWriter&) = delete;
    ParameterWriter(ParameterWriter&&) = delete;
    ParameterWriter& operator=(ParameterWriter&&) = delete;

    auto write() -> void;

   private:
    const mcmc::Result* result_;
    std::unique_ptr<io::detail::TextWriter> writer_;

    auto write_fixed_effects() -> void;
    auto write_random_effects() -> void;
    auto write_residual_variance() -> void;

    auto write_genetic_effect(const GeneticSummary& effect) -> void;

    auto write_summary_statistics(
        std::span<const std::string> terms,
        const PosteriorSummary& stats) -> void;
};

}  // namespace gelex

#endif  // GELEX_IO_MCMC_PARAMETER_WRITER_H_
