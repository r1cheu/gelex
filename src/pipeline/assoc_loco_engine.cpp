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

#include "gelex/pipeline/assoc_loco_engine.h"

#include <Eigen/Core>

#include "detail/assoc_detail.h"
#include "gelex/algo/infer/estimator.h"
#include "gelex/data/genotype/bed_pipe.h"
#include "gelex/data/grm/loco_grm_reader.h"
#include "gelex/data/reader/bim_reader.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/gwas_writer.h"
#include "gelex/model/freq/model.h"
#include "gelex/pipeline/grm_pipe.h"
#include "gelex/pipeline/pheno_pipe.h"
#include "gelex/types/chr_group.h"

namespace gelex
{

AssocLocoEngine::AssocLocoEngine(Config config) : config_(std::move(config)) {}

auto AssocLocoEngine::run(
    PhenoPipe& pheno,
    GrmPipe& grm,
    const df::Index<std::string>& sample_index,
    const AssocObserver& observer,
    const RemlObserver& /*reml_observer*/) -> void
{
    BedPipe bed_pipe(config_.bed_path, sample_index);
    auto bim_path = config_.bed_path;
    auto snp_index
        = std::move(detail::BimReader(bim_path.replace_extension(".bim")))
              .take_info();

    const auto grm_paths = grm.grm_paths();

    FreqModel model(
        std::move(pheno).take_phenotype(),
        std::move(pheno).take_fixed_effects(),
        std::move(grm).take_grms());
    FreqState state(model);

    if (model.genetic().size() != grm_paths.size())
    {
        throw GelexException(
            "Number of genetic components in model does not match number "
            "of GRMs provided.");
    }
    std::vector<LocoGRMReader> loco_readers;
    loco_readers.reserve(grm_paths.size());
    for (const auto& path : grm_paths)
    {
        loco_readers.emplace_back(path, sample_index);
    }

    auto chr_groups = build_chr_groups(true, snp_index);

    notify(
        observer,
        AssocScanSummaryEvent{
            .total_snps = snp_index.size(),
            .chunk_size = config_.chunk_size,
            .loco = true});

    gwas::GwasWriter writer(config_.out_prefix);
    writer.write_header();

    detail::ChrScanner scanner(
        {config_.chunk_size, snp_index.size()}, bed_pipe, observer);

    auto processor = [&](Eigen::Ref<Eigen::MatrixXd> geno, Eigen::VectorXd* f)
    {
        detail::dispatch_assoc_chunk_by_method(
            config_.method, config_.model_type, geno, f);
    };
    auto result_writer = [&](size_t idx, const detail::ChrScanner::SnpResult& r)
    {
        writer.write_result(
            snp_index[idx],
            {.freq = r.freq,
             .beta = r.beta,
             .se = r.se,
             .p_value = r.p_value,
             .pve = r.pve});
    };

    std::vector<LocoRemlResult> loco_results;

    for (const auto& group : chr_groups)
    {
        for (size_t i = 0; i < loco_readers.size(); ++i)
        {
            const auto chr_grm_prefix
                = grm_paths[i].string() + ".chr" + group.name;
            loco_readers[i].load_loco_grm(
                chr_grm_prefix, sample_index, model.genetic()[i].K);
        }

        notify(
            observer,
            AssocLocoPhaseEvent{.chr_name = group.name, .phase = "REML"});

        Estimator estimator(config_.max_iter, config_.tol);
        auto v_inv = estimator.fit(model, state);
        detail::update_assoc_input(
            scanner.assoc_input(), model, state, std::move(v_inv));

        LocoRemlResult result;
        result.chr_name = group.name;
        result.loglike = estimator.loglike();
        result.converged = estimator.is_converged();
        result.residual_variance = state.residual().variance;
        for (const auto& g : state.genetic())
        {
            result.genetic.push_back(
                {.type = g.type,
                 .variance = g.variance,
                 .heritability = g.heritability});
        }
        loco_results.push_back(std::move(result));

        notify(
            observer,
            AssocLocoPhaseEvent{.chr_name = group.name, .phase = "SCAN"});

        scanner.scan(group, processor, result_writer);
    }

    writer.finalize();

    notify(observer, AssocLocoRemlSummaryEvent{.results = loco_results});

    notify(observer, AssocCompleteEvent{.out_prefix = config_.out_prefix});
}

}  // namespace gelex
