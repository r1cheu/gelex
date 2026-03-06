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

#include "gelex/pipeline/assoc_normal_engine.h"

#include <Eigen/Core>

#include "assoc_detail.h"
#include "gelex/algo/infer/estimator.h"
#include "gelex/data/genotype/bed_pipe.h"
#include "gelex/data/loader/bim_loader.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/freq/model.h"
#include "gelex/pipeline/grm_pipe.h"
#include "gelex/pipeline/pheno_pipe.h"
#include "gelex/pipeline/report/gwas_writer.h"
#include "gelex/types/chr_group.h"

namespace gelex
{

AssocNormalEngine::AssocNormalEngine(Config config) : config_(std::move(config))
{
}

auto AssocNormalEngine::run(
    PhenoPipe& pheno,
    GrmPipe& grm,
    const AssocObserver& observer,
    const RemlObserver& reml_observer) -> void
{
    BedPipe bed_pipe(config_.bed_path, pheno.sample_manager());
    auto bim_path = config_.bed_path;
    auto snp_effects
        = std::move(detail::BimLoader(bim_path.replace_extension(".bim")))
              .take_info();

    FreqModel model(pheno, grm);
    FreqState state(model);
    Estimator estimator(config_.max_iter, config_.tol, reml_observer);

    notify(observer, AssocRemlStartedEvent{.chr_name = ""});

    auto v_inv = estimator.fit(model, state);
    auto chr_groups = build_chr_groups(false, snp_effects);

    notify(
        observer,
        AssocScanSummaryEvent{
            .total_snps = snp_effects.size(),
            .chunk_size = config_.chunk_size,
            .loco = false});

    gwas::GwasWriter writer(config_.out_prefix);
    writer.write_header();

    detail::ChrScanner scanner(
        {config_.chunk_size, snp_effects.size()}, bed_pipe, observer);
    detail::update_assoc_input(
        scanner.assoc_input(), model, state, std::move(v_inv));

    auto processor = [&](Eigen::Ref<Eigen::MatrixXd> geno, Eigen::VectorXd* f)
    {
        detail::dispatch_assoc_chunk_by_method(
            config_.method, config_.model_type, geno, f);
    };
    auto result_writer = [&](size_t idx, const detail::ChrScanner::SnpResult& r)
    {
        writer.write_result(
            snp_effects[idx],
            {.freq = r.freq, .beta = r.beta, .se = r.se, .p_value = r.p_value});
    };
    for (const auto& group : chr_groups)
    {
        scanner.scan(group, processor, result_writer);
    }

    writer.finalize();

    notify(observer, AssocCompleteEvent{.out_prefix = config_.out_prefix});
}

}  // namespace gelex
