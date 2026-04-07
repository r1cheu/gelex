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

#ifndef GELEX_POST_GENETIC_VARIANCE_PROCESSOR_H_
#define GELEX_POST_GENETIC_VARIANCE_PROCESSOR_H_

#include <span>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/stats/diagnostics.h"
#include "gelex/infra/logging/post_event.h"
#include "gelex/io/binary_reader.h"
#include "gelex/model/bayes/genotype_storage.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

struct GeneticInput
{
    const bayes::GenotypeStorage* genotype;
    GeneticKind kind;
};

class GeneticVarianceProcessor
{
   public:
    GeneticVarianceProcessor(
        const GeneticInput& input,
        std::span<const detail::BinaryReader> readers);

    auto process(
        std::size_t chain_idx,
        Eigen::Index col_begin,
        Eigen::Index col_end) -> const Eigen::MatrixXd&;

    [[nodiscard]] auto last_variances() const
        -> Eigen::Ref<const Eigen::RowVectorXd>
    {
        return last_variances_;
    }

    [[nodiscard]] auto has_components() const -> bool
    {
        return n_components_ > 0;
    }

    [[nodiscard]] auto build_diagnostics(double hdpi_threshold) const
        -> std::vector<ParameterDiag>;

   private:
    Eigen::Ref<const Eigen::MatrixXd> matrix_;
    GeneticKind kind_;
    Eigen::Index n_components_;
    std::vector<Eigen::Map<const Eigen::MatrixXd>> coeff_maps_;
    std::vector<Eigen::Map<const Eigen::MatrixX<int8_t>>> tracker_maps_;
    Chains component_variance_chains_;
    Eigen::MatrixXd gebv_chunk_;
    Eigen::MatrixXd component_gebv_chunk_;
    Eigen::MatrixXd masked_beta_chunk_;
    Eigen::RowVectorXd last_variances_;
};

}  // namespace gelex

#endif  // GELEX_POST_GENETIC_VARIANCE_PROCESSOR_H_
