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

#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <fmt/format.h>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "gelex/algo/gwas/assoc_tester.h"
#include "gelex/algo/gwas/assoc_type.h"
#include "gelex/algo/reml/estimator.h"
#include "gelex/algo/reml/operators.h"
#include "gelex/algo/reml/summary.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/marker_range.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/freq/design.h"
#include "gelex/freq/model.h"
#include "gelex/io/gwas_writer.h"
#include "gelex/io/loco_reader.h"
#include "gelex/types/fixed_designs.h"

#include "cli/cli_helper.h"
#include "cli/common_data.h"
#include "cli/formatter.h"
#include "cli/reml_data.h"
#include "cli/reml_reporter.h"
#include "cli/report_printer.h"
#include "reporter.h"

class AssocDataHandler
{
   public:
    AssocDataHandler(const cli::AssocConfig& config, gelex::Bed& bed) noexcept
        : bed_(bed), loader_(config.random)
    {
    }

    auto load_indices(
        std::vector<const gelex::DataFrameIndex<std::string>*>& indices) -> void
    {
        indices.push_back(&bed_.sample_index());
        loader_.load_indices(indices);
        for (const auto& grm_index : loader_.grm_indices())
        {
            cli::printer().line("   GRM        : {} samples", grm_index.size());
        }
    }

    auto gather(const gelex::DataFrameIndex<std::string>& common_index) -> void
    {
        bed_.gather(common_index);
        loader_.gather(common_index);
    }

    auto results() && -> std::vector<gelex::freq::RandomDesign>
    {
        return std::move(loader_).results();
    }

   private:
    gelex::Bed& bed_;
    cli::RemlDataLoader loader_;
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

    cli::printer().block(gelex::section("Dataset Summary:"));

    auto bed = gelex::open_bed(config.bfile);
    AssocDataHandler handler(config, bed);
    cli::BaseData data = cli::load_base_data(handler, config.base_data);
    auto random_designs = std::move(handler).results();

    const auto& sample_index = bed.sample_index();
    cli::printer().line(
        "   Intersection : {} common samples", sample_index.size());
    if (sample_index.size() == 0)
    {
        throw gelex::GelexException(
            "No common samples across phenotype, genotype (.fam), GRM, and "
            "covariates. Check that sample IDs match across input files.");
    }

    const auto& bim = bed.bim();

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
        = [&](const gelex::MarkerRange& range, const gelex::GwasOperators& ops)
    {
        const auto n_samples = ops.n_samples();

        for (auto start = range.start; start < range.end;
             start += static_cast<Eigen::Index>(config.chunk_size))
        {
            const auto end = std::min(
                start + static_cast<Eigen::Index>(config.chunk_size),
                range.end);
            const auto current_chunk_size = end - start;

            tester->resize(n_samples, current_chunk_size);
            bed.read_into<double>(tester->genotype_buffer(), start);

            auto results = tester->run(ops);
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

        auto fit = estimator.fit(model, state);
        reml_reporter.show_result(model, fit.summary, config.max_iter);

        reporter.start_scan(total_snps, config.chunk_size, false);

        scan_range(
            {"all", 0, static_cast<Eigen::Index>(bim.rows())}, fit.operators);
    }
    else
    {
        std::vector<Eigen::MatrixXd> whole_grms;
        whole_grms.reserve(config.random.grm.size());
        std::vector<gelex::LocoReader> loco_readers;
        loco_readers.reserve(config.random.grm.size());
        for (const auto& path : config.random.grm)
        {
            whole_grms.push_back(gelex::read_grm(path, &sample_index, false));
            loco_readers.emplace_back(whole_grms.back());
        }

        std::vector<Eigen::Index> grm_slots;
        for (auto&& [idx, design] : std::views::enumerate(model.random()))
        {
            if (design.kind == gelex::freq::RandomKind::Grm)
            {
                grm_slots.push_back(static_cast<Eigen::Index>(idx));
            }
        }

        auto ranges = gelex::chromosome_ranges(bim);

        reporter.start_scan(total_snps, config.chunk_size, true);

        std::vector<gelex::LocoRemlResult> loco_results;

        gelex::Estimator estimator(config.max_iter, config.tolerance);
        for (const auto& range : ranges)
        {
            for (std::size_t i = 0; i < loco_readers.size(); ++i)
            {
                const auto chr_grm_prefix = fmt::format(
                    "{}.chr{:02d}",
                    config.random.grm[i],
                    std::stoi(range.label));
                loco_readers[i].load_into(
                    chr_grm_prefix,
                    sample_index,
                    model.random()[grm_slots[i]].K);
            }

            reporter.show_loco_phase(range.label, "REML");

            auto fit = estimator.fit(model, state);

            loco_results.push_back({range.label, std::move(fit.summary)});

            reporter.show_loco_phase(range.label, "SCAN");

            scan_range(range, fit.operators);
        }

        reporter.show_loco_reml_summary(loco_results);
    }

    reporter.show_complete(config.out);

    return 0;
}
