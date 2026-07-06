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

#include "command.h"

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "cli/cli_helper.h"
#include "cli/common_data.h"
#include "cli/formatter.h"
#include "cli/reml_reporter.h"
#include "cli/report_printer.h"
#include "gelex/algo/gwas/assoc_tester.h"
#include "gelex/algo/gwas/assoc_type.h"
#include "gelex/algo/reml/estimator.h"
#include "gelex/algo/reml/loco_result.h"
#include "gelex/algo/reml/result.h"
#include "gelex/data/bed.h"
#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/marker_range.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/freq/model.h"
#include "gelex/io/gwas/writer.h"
#include "gelex/io/loco_reader.h"
#include "gelex/types/fixed_designs.h"
#include "reporter.h"

class AssocDataHandler
{
   public:
    explicit AssocDataHandler(const cli::AssocConfig& config) noexcept
        : config_(config)
    {
    }

    auto load_indices(std::vector<gelex::DataFrameIndex<std::string>*>& indices)
        -> void
    {
        fam_index_ = gelex::read_fam(config_.bfile + ".fam").index();
        indices.push_back(&fam_index_);

        grm_indices_.reserve(config_.grm.size());
        for (const auto& path : config_.grm)
        {
            grm_indices_.emplace_back(gelex::read_grm_ids(path));
            cli::printer().line(
                "   GRM        : {} samples", grm_indices_.back().size());
            indices.push_back(&grm_indices_.back());
        }
    }

    auto gather(const gelex::DataFrameIndex<std::string>& common_index) -> void
    {
        sample_index_ = common_index;
        random_designs_ = gelex::make_grm_designs(config_.grm, common_index);
    }

    auto results() && -> std::pair<
        gelex::DataFrameIndex<std::string>,
        std::vector<gelex::freq::RandomDesign>>
    {
        return {std::move(sample_index_), std::move(random_designs_)};
    }

   private:
    const cli::AssocConfig& config_;
    gelex::DataFrameIndex<std::string> fam_index_;
    gelex::DataFrameIndex<std::string> sample_index_;
    std::vector<gelex::DataFrameIndex<std::string>> grm_indices_;
    std::vector<gelex::freq::RandomDesign> random_designs_;
};

auto assoc_execute(const cli::AssocConfig& config) -> int
{
    cli::setup_parallelization(config.threads);

    const bool is_joint = config.mode.size() > 1;
    const gelex::AssocType test_type
        = is_joint ? gelex::AssocType::Joint : gelex::AssocType::Single;
    const gelex::GeneticMode mode
        = is_joint ? gelex::GeneticMode::A : *config.mode.each().begin();
    const gelex::GenotypeMethod geno_method{config.geno_method};

    cli::AssocReporter reporter;

    cli::printer().block(gelex::section("[Dataset Summary]"));

    AssocDataHandler handler(config);
    cli::BaseData data = cli::load_base_data(handler, config.base_data);
    auto [sample_index, random_designs] = std::move(handler).results();

    cli::printer().line(
        "   Intersection : {} common samples", sample_index.size());
    if (sample_index.size() == 0)
    {
        throw gelex::GelexException(
            "No common samples across phenotype, genotype (.fam), GRM, and "
            "covariates. Check that sample IDs match across input files.");
    }

    auto bed = gelex::open_bed(config.bfile, sample_index);
    auto bim = gelex::read_bim(config.bfile + ".bim");

    gelex::FreqModel model(
        std::move(data.phenotype),
        std::move(data.fixed_design),
        std::move(random_designs));
    gelex::FreqState state(model);

    auto tester = gelex::AssocTester::make(test_type, mode, geno_method);
    gelex::GwasWriter writer(config.out, bim, test_type);

    const auto total_snps = static_cast<std::size_t>(bim.rows());
    std::size_t progress = 0;

    const auto scan_range
        = [&](const gelex::MarkerRange& range, const gelex::RemlResult& reml)
    {
        const auto n_samples = reml.n_samples();

        for (auto start = range.start; start < range.end;
             start += static_cast<Eigen::Index>(config.chunk_size))
        {
            const auto end = std::min(
                start + static_cast<Eigen::Index>(config.chunk_size),
                range.end);
            const auto current_chunk_size = end - start;

            tester->resize(n_samples, current_chunk_size);
            bed.read_into<double>(tester->genotype_buffer(), start);

            auto results = tester->run(reml);
            writer.write(static_cast<std::size_t>(start), results);

            progress += static_cast<std::size_t>(current_chunk_size);
            reporter.update_scan_progress(progress, total_snps);
        }
    };

    cli::RemlReporter reml_reporter;

    if (!config.loco)
    {
        gelex::Estimator estimator(
            config.max_iter, config.tolerance, reml_reporter.as_observer());

        reporter.show_reml_started("");

        auto reml = estimator.fit(model, state);
        reml_reporter.show_result(
            model,
            state,
            estimator.is_converged(),
            estimator.iter_count(),
            config.max_iter,
            estimator.loglike());

        reporter.start_scan(total_snps, config.chunk_size, false);

        scan_range({"all", 0, static_cast<Eigen::Index>(bim.rows())}, reml);
    }
    else
    {
        std::vector<gelex::LocoReader> loco_readers;
        loco_readers.reserve(config.grm.size());
        for (const auto& path : config.grm)
        {
            loco_readers.emplace_back(path, sample_index);
        }

        auto ranges = gelex::chromosome_ranges(bim);

        reporter.start_scan(total_snps, config.chunk_size, true);

        std::vector<gelex::LocoRemlResult> loco_results;

        for (const auto& range : ranges)
        {
            for (std::size_t i = 0; i < loco_readers.size(); ++i)
            {
                const auto chr_grm_prefix = fmt::format(
                    "{}.chr{:02d}", config.grm[i], std::stoi(range.label));
                loco_readers[i].load_into(
                    chr_grm_prefix, sample_index, model.random()[i].K);
            }

            reporter.show_loco_phase(range.label, "REML");

            gelex::Estimator estimator(config.max_iter, config.tolerance);
            auto reml = estimator.fit(model, state);

            {
                gelex::LocoRemlResult r;
                r.chr_name = range.label;
                r.loglike = estimator.loglike();
                r.converged = estimator.is_converged();
                r.residual_variance = state.residual().variance;
                for (std::size_t i = 0; i < state.random().size(); ++i)
                {
                    const auto& random = state.random()[i];
                    r.random.push_back(
                        {.name = model.random()[i].name,
                         .variance = random.variance,
                         .variance_ratio = random.variance_ratio});
                }
                loco_results.push_back(std::move(r));
            }

            reporter.show_loco_phase(range.label, "SCAN");

            scan_range(range, reml);
        }

        reporter.show_loco_reml_summary(loco_results);
    }

    reporter.show_complete(config.out);

    return 0;
}
