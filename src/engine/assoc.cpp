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

#include "gelex/engine/assoc.h"

#include <algorithm>
#include <cstddef>

#include <Eigen/Core>
#include <string>
#include <utility>
#include <vector>

#include "gelex/algo/gwas/assoc_tester.h"
#include "gelex/algo/reml/estimator.h"
#include "gelex/algo/reml/result.h"
#include "gelex/data/chr_group.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype/bed_pipe.h"
#include "gelex/data/pipe/grm.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/assoc_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/logging/reml_event.h"
#include "gelex/io/grm/loco_reader.h"
#include "gelex/io/gwas/writer.h"
#include "gelex/model/freq/model.h"

namespace gelex
{

AssocEngine::AssocEngine(Config config) : config_(std::move(config)) {}

auto AssocEngine::run(
    PhenoPipe& pheno,
    GrmPipe& grm,
    const dataframe::Index<std::string>& sample_index,
    const AssocObserver& observer,
    const RemlObserver& reml_observer) -> void
{
    genotype::BedPipe bed_pipe(config_.bfile_prefix, sample_index);
    auto bim = read_bim(config_.bfile_prefix + ".bim");

    auto phenotype = std::move(pheno).take_phenotype();
    auto fixed_design = std::move(pheno).take_fixed_design();
    auto grm_paths = grm.grm_paths();
    auto grms = std::move(grm).take_grms();

    FreqModel model(
        std::move(phenotype), std::move(fixed_design), std::move(grms));
    FreqState state(model);

    auto tester
        = AssocTester::make(config_.test_type, config_.mode, config_.method);
    gwas::GwasWriter writer(config_.out_prefix, bim, config_.test_type);

    const auto total_snps = static_cast<size_t>(bim.rows());
    size_t progress = 0;

    const auto scan_groups
        = [&](const std::vector<ChrGroup>& groups, const RemlResult& reml)
    {
        const auto n_samples = reml.n_samples();
        const auto eigen_chunk = static_cast<Eigen::Index>(config_.chunk_size);

        for (const auto& group : groups)
        {
            for (const auto& [range_start, range_end] : group.ranges)
            {
                for (auto start = range_start; start < range_end;
                     start += eigen_chunk)
                {
                    const auto end = std::min(start + eigen_chunk, range_end);
                    const auto current_chunk_size = end - start;

                    tester->resize(n_samples, current_chunk_size);
                    bed_pipe.load_chunk(tester->genotype_buffer(), start, end);

                    auto results = tester->run(reml);
                    writer.write(static_cast<size_t>(start), results);

                    progress += static_cast<size_t>(current_chunk_size);
                    notify(
                        observer,
                        AssocScanProgressEvent{
                            .current = progress, .total = total_snps});
                }
            }
        }
    };

    if (!config_.loco)
    {
        reml::Estimator estimator(config_.max_iter, config_.tol, reml_observer);

        notify(observer, AssocRemlStartedEvent{.chr_name = ""});

        auto reml = estimator.fit(model, state);

        auto chr_groups = build_chr_groups(false, bim);

        notify(
            observer,
            AssocScanSummaryEvent{
                .total_snps = total_snps,
                .chunk_size = config_.chunk_size,
                .loco = false});

        scan_groups(chr_groups, reml);
    }
    else
    {
        if (model.genetic().size() != grm_paths.size())
        {
            throw GelexException(
                "Number of genetic components in model does not match "
                "number of GRMs provided.");
        }

        std::vector<LocoReader> loco_readers;
        loco_readers.reserve(grm_paths.size());
        for (const auto& path : grm_paths)
        {
            loco_readers.emplace_back(path, sample_index);
        }

        auto chr_groups = build_chr_groups(true, bim);

        notify(
            observer,
            AssocScanSummaryEvent{
                .total_snps = total_snps,
                .chunk_size = config_.chunk_size,
                .loco = true});

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

            reml::Estimator estimator(config_.max_iter, config_.tol);
            auto reml = estimator.fit(model, state);

            {
                LocoRemlResult r;
                r.chr_name = group.name;
                r.loglike = estimator.loglike();
                r.converged = estimator.is_converged();
                r.residual_variance = state.residual().variance;
                for (const auto& g : state.genetic())
                {
                    r.genetic.push_back(
                        {.type = g.type,
                         .variance = g.variance,
                         .heritability = g.heritability});
                }
                loco_results.push_back(std::move(r));
            }

            notify(
                observer,
                AssocLocoPhaseEvent{.chr_name = group.name, .phase = "SCAN"});

            scan_groups({group}, reml);
        }

        notify(observer, AssocLocoRemlSummaryEvent{.results = loco_results});
    }

    notify(observer, AssocCompleteEvent{.out_prefix = config_.out_prefix});
}

}  // namespace gelex
