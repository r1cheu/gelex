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

#ifndef GELEX_POST_HERITABILITY_H_
#define GELEX_POST_HERITABILITY_H_

#include <span>
#include <vector>

#include "gelex/infra/logging/post_event.h"
#include "gelex/infra/stats/diagnostics.h"
#include "gelex/io/binary_reader.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class HeritabilityPosteriorProcessor
{
   public:
    explicit HeritabilityPosteriorProcessor(
        std::span<const io::BinaryReader> readers,
        double hdpi_threshold,
        const stats::Chains& genetic_variances,
        std::span<const GeneticMode> kinds);

    [[nodiscard]] auto process() -> std::vector<ParameterDiag>;

   private:
    auto assemble_phenotypic_variance() const -> stats::Chains;

    std::span<const io::BinaryReader> readers_;
    double hdpi_threshold_;
    const stats::Chains& genetic_variances_;
    std::span<const GeneticMode> kinds_;
};

}  // namespace gelex

#endif  // GELEX_POST_HERITABILITY_H_
