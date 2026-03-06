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

#include "assoc_detail.h"

#include "gelex/algo/gwas/association_test.h"
#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/freq/model.h"

namespace gelex::detail
{

auto dispatch_assoc_chunk_by_method(
    GenotypeProcessMethod method,
    ModelType model_type,
    Eigen::Ref<Eigen::MatrixXd> genotype,
    Eigen::VectorXd* freqs) -> void
{
    if (!is_center_family_method(method))
    {
        throw InvalidInputException(
            "assoc --geno-method supports only center-family methods: "
            "2 (center-hwe), 4 (orth-center-hwe), 6 (center), 8 (orth-center)");
    }
    if (model_type == ModelType::A)
    {
        process_matrix<GeneticEffectType::Add>(method, genotype, freqs);
    }
    else
    {
        process_matrix<GeneticEffectType::Dom>(method, genotype, freqs);
    }
}

auto update_assoc_input(
    AssocInput& input,
    const FreqModel& model,
    const FreqState& state,
    Eigen::MatrixXd&& v_inv) -> void
{
    input.V_inv = std::move(v_inv);
    input.V_inv_y
        = input.V_inv
          * (model.phenotype() - model.fixed().X * state.fixed().coeff);
}

ChrScanner::ChrScanner(Config config, BedPipe& bed, AssocObserver observer)
    : config_(config), bed_(bed), observer_(std::move(observer))
{
    const auto cs = static_cast<Eigen::Index>(config_.chunk_size);
    output_.resize(cs);
    freqs_.resize(cs);
}

auto ChrScanner::scan(
    const ChrGroup& group,
    const GenoProcessor& processor,
    const ResultWriter& writer) -> void
{
    const auto n_samples = input_.V_inv.rows();
    const auto cs = static_cast<size_t>(config_.chunk_size);

    for (const auto& [range_start, range_end] : group.ranges)
    {
        for (auto start = range_start; start < range_end;
             start += static_cast<Eigen::Index>(cs))
        {
            const auto end
                = std::min(start + static_cast<Eigen::Index>(cs), range_end);
            const auto current_chunk_size = end - start;

            input_.resize(n_samples, current_chunk_size);
            output_.resize(current_chunk_size);
            freqs_.resize(current_chunk_size);

            bed_.load_chunk(input_.Z, start, end);

            processor(input_.Z, &freqs_);
            gwas::wald_test(input_, output_);

            for (Eigen::Index i = 0; i < current_chunk_size; ++i)
            {
                writer(
                    static_cast<size_t>(start + i),
                    {.freq = freqs_(i),
                     .beta = output_.beta(i),
                     .se = output_.se(i),
                     .p_value = output_.p_value(i)});
            }

            progress_counter_ += static_cast<size_t>(current_chunk_size);

            notify(
                observer_,
                AssocScanProgressEvent{
                    .current = progress_counter_, .total = config_.total_snps});
        }
    }
}

}  // namespace gelex::detail
