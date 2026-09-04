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
#include <fmt/ranges.h>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/encode/encoder.h"
#include "gelex/data/fixed_design.h"
#include "gelex/data/grm/io.h"
#include "gelex/data/marker_range.h"
#include "gelex/exception.h"
#include "gelex/freq/design.h"
#include "gelex/freq/gwas/assoc_tester.h"
#include "gelex/freq/gwas/assoc_type.h"
#include "gelex/freq/gwas/loco_grm.h"
#include "gelex/freq/gwas/writer.h"
#include "gelex/freq/model.h"
#include "gelex/freq/reml/estimator.h"
#include "gelex/freq/reml/operators.h"
#include "gelex/freq/reml/summary.h"

#include "cli/assoc/progress.h"
#include "cli/common_data.h"
#include "cli/formatter.h"
#include "cli/reml_data.h"
#include "cli/reml_reporter.h"
#include "cli/report_printer.h"
#include "cli/runtime.h"
#include "cli/summary.h"

namespace
{

auto make_scan_summary(std::size_t total_snps, int chunk_size, bool loco)
    -> cli::Summary
{
    cli::Summary summary{"Association Scan"};
    summary.field("Variants", "{}", total_snps)
        .field("Chunk size", "{}", chunk_size);
    if (loco)
    {
        summary.field("Mode", "LOCO");
    }
    return summary;
}

}  // namespace

class AssocDataHandler
{
   public:
    AssocDataHandler(const cli::AssocConfig& config, gelex::Bed& bed)
        : bed_(&bed), loader_(config.random)
    {
    }

    auto load_indices(
        std::vector<const gelex::DataFrameIndex<std::string>*>& indices) -> void
    {
        indices.push_back(&bed_->sample_index());
        loader_.load_indices(indices);
    }

    auto gather(const gelex::DataFrameIndex<std::string>& common_index) -> void
    {
        bed_->gather(common_index);
        loader_.gather(common_index);
    }

    auto results() && -> std::vector<gelex::freq::RandomDesign>
    {
        return std::move(loader_).results();
    }

   private:
    gelex::Bed* bed_;
    cli::RemlDataLoader loader_;
};

auto assoc_execute(const cli::AssocConfig& config) -> int
{
    if (!config.random.has_random_effect())
    {
        throw gelex::GelexException(
            "association testing needs at least one random effect; provide "
            "--grm, --drand, --qrand, or --interaction");
    }
    if (config.loco && !config.random.interactions.empty())
    {
        throw gelex::GelexException(
            "--loco does not support --interaction: an interaction kernel "
            "stays whole-genome and cannot be left-one-chromosome-out");
    }

    cli::setup_parallelization(config.threads);

    const bool is_joint = config.mode.size() > 1;
    const gelex::AssocType test_type
        = is_joint ? gelex::AssocType::Joint : gelex::AssocType::Single;
    const gelex::GeneticMode mode
        = is_joint ? gelex::GeneticMode::A : *config.mode.each().begin();
    const gelex::GenotypeMethod geno_method{config.geno_method};

    auto bed = gelex::open_bed(config.bfile);
    AssocDataHandler handler(config, bed);
    cli::BaseData data = cli::load_base_data(handler, config.base_data);
    auto random_designs = std::move(handler).results();

    const auto& sample_index = bed.sample_index();
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

    cli::Summary{"Dataset Summary"}
        .field("Trait", "{}", data.pheno_name)
        .field("Samples", "{}", model.num_individuals())
        .field("Variants", "{}", bed.num_snps())
        .show();

    const auto random_effect_names
        = model.random()
          | std::views::transform([](const auto& design)
                                  { return design.name; });
    cli::Summary{"Model Summary"}
        .field("Fixed terms", "{}", model.fixed().column_names().size())
        .field("Random effects", "{}", fmt::join(random_effect_names, ", "))
        .show();

    gelex::FreqState state(model);

    auto tester = gelex::AssocTester::make(test_type, mode, geno_method);
    const gelex::LocusEncoder encoder{bed};
    gelex::GwasWriter writer(config.out, bim, test_type);

    const auto total_snps = static_cast<std::size_t>(bim.rows());
    std::size_t scanned_snps = 0;

    const auto scan_range = [&](const gelex::MarkerRange& range,
                                const gelex::GwasOperators& ops,
                                cli::AssocProgress& progress)
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
            auto results = tester->run(encoder, start, ops);
            writer.write(static_cast<std::size_t>(start), results);

            scanned_snps += static_cast<std::size_t>(current_chunk_size);
            progress(scanned_snps);
        }
    };

    cli::RemlReporter reml_reporter;

    if (!config.loco)
    {
        gelex::Estimator estimator(
            config.max_iter, config.tolerance, reml_reporter.as_observer());

        auto fit = estimator.fit(model, state);
        reml_reporter.show_result(model, fit.summary, config.max_iter);

        make_scan_summary(total_snps, config.chunk_size, false).show();
        cli::AssocProgress progress{total_snps};

        scan_range(
            {"all", 0, static_cast<Eigen::Index>(bim.rows())},
            fit.operators,
            progress);
        progress.finish();
    }
    else
    {
        std::vector<Eigen::MatrixXd> whole_grms;
        whole_grms.reserve(config.random.grm.size());
        std::vector<gelex::LocoGrmBuilder> loco_builders;
        loco_builders.reserve(config.random.grm.size());
        for (const auto& path : config.random.grm)
        {
            whole_grms.push_back(gelex::read_grm(path, &sample_index, false));
            loco_builders.emplace_back(whole_grms.back());
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

        make_scan_summary(total_snps, config.chunk_size, true).show();
        cli::AssocProgress progress{total_snps};

        std::vector<gelex::LocoRemlResult> loco_results;

        gelex::Estimator estimator(config.max_iter, config.tolerance);
        for (const auto& range : ranges)
        {
            progress.start_reml();

            for (std::size_t i = 0; i < loco_builders.size(); ++i)
            {
                const auto chr_grm_prefix = fmt::format(
                    "{}.chr{:02d}",
                    config.random.grm[i],
                    std::stoi(range.label));
                const auto chromosome_grm
                    = gelex::read_grm(chr_grm_prefix, &sample_index, false);
                loco_builders[i].build_into(
                    chromosome_grm, model.random()[grm_slots[i]].K);
            }

            auto fit = estimator.fit(model, state);

            loco_results.push_back({range.label, std::move(fit.summary)});

            progress.start_scan();

            scan_range(range, fit.operators, progress);
        }

        progress.finish();
        cli::print_loco_reml_summary(loco_results);
    }

    cli::printer().block(cli::results_saved(config.out, ".gwas.tsv, .log"));

    return 0;
}
