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
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>
#include <Eigen/Core>

#include "cli/cli_helper.h"
#include "cli/common_data.h"
#include "cli/reml_reporter.h"
#include "cli/report_printer.h"
#include "gelex/algo/gwas/assoc_tester.h"
#include "gelex/algo/gwas/assoc_type.h"
#include "gelex/algo/reml/estimator.h"
#include "gelex/algo/reml/loco_result.h"
#include "gelex/algo/reml/result.h"
#include "gelex/data/bed.h"
#include "gelex/data/chr_group.h"
#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/freq/model.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/logging/formatter.h"
#include "gelex/infra/stats/rank_inverse_norm_transform.h"
#include "gelex/io/grm/loco_reader.h"
#include "gelex/io/gwas/joint_cov_writer.h"
#include "gelex/io/gwas/writer.h"
#include "gelex/types/fixed_designs.h"
#include "reporter.h"

struct AssocData
{
    gelex::dataframe::Index<std::string> sample_index;
    std::vector<std::string> grm_prefixes;
    std::vector<gelex::freq::RandomDesign> random_designs;
};

class AssocDataHandler
{
   public:
    auto load_indices(
        CLI::App& cmd,
        std::vector<gelex::dataframe::Index<std::string>*>& indices) -> void
    {
        fam_index_ = gelex::read_fam(
                         cmd.get_option("--bfile")->as<std::string>() + ".fam")
                         .index();
        indices.push_back(&fam_index_);

        grm_prefixes_ = cmd.get_option("--grm")->as<std::vector<std::string>>();
        grm_indices_.reserve(grm_prefixes_.size());
        for (const auto& path : grm_prefixes_)
        {
            grm_indices_.emplace_back(gelex::read_grm_ids(path));
            cli::printer().line(
                "   GRM        : {} samples", grm_indices_.back().size());
            indices.push_back(&grm_indices_.back());
        }
    }

    auto gather(const gelex::dataframe::Index<std::string>& common_index)
        -> void
    {
        sample_index_ = common_index;
        random_designs_ = gelex::make_grm_designs(grm_prefixes_, common_index);
    }

    auto results() && -> AssocData
    {
        return AssocData{
            .sample_index = std::move(sample_index_),
            .grm_prefixes = std::move(grm_prefixes_),
            .random_designs = std::move(random_designs_)};
    }

   private:
    gelex::dataframe::Index<std::string> fam_index_;
    gelex::dataframe::Index<std::string> sample_index_;
    std::vector<std::string> grm_prefixes_;
    std::vector<gelex::dataframe::Index<std::string>> grm_indices_;
    std::vector<gelex::freq::RandomDesign> random_designs_;
};

auto assoc_execute(CLI::App& cmd) -> int
{
    cli::setup_parallelization(cmd.get_option("--threads")->as<int>());

    const gelex::AssocType test_type
        = cmd.get_option("--test")->as<std::string>() == "joint"
              ? gelex::AssocType::Joint
              : gelex::AssocType::Single;
    const gelex::GeneticMode mode
        = test_type == gelex::AssocType::Single
                  && cmd.get_option("--mode")->as<std::string>() == "D"
              ? gelex::GeneticMode::D
              : gelex::GeneticMode::A;
    const gelex::GenotypeMethod geno_method{cli::parse_genotype_method(
        cmd.get_option("--geno-method")->as<std::string>())};
    const bool write_cov{cmd.get_option("--write-cov")->count() > 0};

    if (write_cov
        && (test_type != gelex::AssocType::Joint
            || geno_method != gelex::GenotypeMethod::Center))
    {
        throw gelex::GelexException(
            "assoc --write-cov requires --test joint --geno-method C");
    }

    cli::AssocReporter reporter;

    reporter.show_banner();
    reporter.show_config(
        mode,
        test_type,
        cmd.get_option("--loco")->count() > 0,
        geno_method,
        cmd.get_option("--max-iter")->as<int>(),
        cmd.get_option("--tol")->as<double>());

    cli::printer().block(gelex::section("[Dataset Summary]"));

    AssocDataHandler handler;
    cli::BaseData data = cli::load_base_data(handler, cmd);
    auto assoc_data = std::move(handler).results();

    cli::printer().line(
        "   Intersection : {} common samples", assoc_data.sample_index.size());
    if (assoc_data.sample_index.size() == 0)
    {
        throw gelex::GelexException(
            "No common samples across phenotype, genotype (.fam), GRM, and "
            "covariates. Check that sample IDs match across input files.");
    }

    if (cmd.get_option("--transform")->as<std::string>() != "none")
    {
        gelex::stats::RankInverseNormTransform transformer(
            cmd.get_option("--int-offset")->as<double>());
        auto logger = gelex::logging::get();

        if (cmd.get_option("--transform")->as<std::string>() == "dint")
        {
            logger->info(
                "   Method: Direct INT (DINT), offset (k): {}",
                cmd.get_option("--int-offset")->as<double>());
            transformer.apply_dint(data.phenotype);
        }
        else if (cmd.get_option("--transform")->as<std::string>() == "iint")
        {
            logger->info(
                "   Method: Indirect INT (IINT), offset (k): {}",
                cmd.get_option("--int-offset")->as<double>());
            transformer.apply_iint(data.phenotype, data.fixed_design.X);
            data.fixed_design = gelex::FixedDesign::make(data.phenotype.size());
        }
    }

    auto bed = gelex::open_bed(
        cmd.get_option("--bfile")->as<std::string>(), assoc_data.sample_index);
    auto bim = gelex::read_bim(
        cmd.get_option("--bfile")->as<std::string>() + ".bim");

    gelex::FreqModel model(
        std::move(data.phenotype),
        std::move(data.fixed_design),
        std::move(assoc_data.random_designs));
    gelex::FreqState state(model);

    auto tester = gelex::AssocTester::make(test_type, mode, geno_method);
    gelex::gwas::GwasWriter writer(
        cmd.get_option("--out")->as<std::string>(), bim, test_type);
    std::optional<gelex::gwas::JointCovWriter> joint_cov_writer;
    if (write_cov)
    {
        joint_cov_writer.emplace(
            cmd.get_option("--out")->as<std::string>(), bim);
    }

    const auto total_snps = static_cast<std::size_t>(bim.rows());
    std::size_t progress = 0;

    const auto scan_groups = [&](const std::vector<gelex::ChrGroup>& groups,
                                 const gelex::RemlResult& reml)
    {
        const auto n_samples = reml.n_samples();

        for (const auto& group : groups)
        {
            for (const auto& [range_start, range_end] : group.ranges)
            {
                for (auto start = range_start; start < range_end;
                     start += static_cast<Eigen::Index>(
                         cmd.get_option("--chunk-size")->as<int>()))
                {
                    const auto end = std::min(
                        start
                            + static_cast<Eigen::Index>(
                                cmd.get_option("--chunk-size")->as<int>()),
                        range_end);
                    const auto current_chunk_size = end - start;

                    tester->resize(n_samples, current_chunk_size);
                    bed.read_into<double>(tester->genotype_buffer(), start);

                    auto results = tester->run(reml);
                    writer.write(static_cast<std::size_t>(start), results);
                    if (joint_cov_writer)
                    {
                        joint_cov_writer->write(
                            static_cast<std::size_t>(start), results);
                    }

                    progress += static_cast<std::size_t>(current_chunk_size);
                    reporter.update_scan_progress(progress, total_snps);
                }
            }
        }
    };

    cli::RemlReporter reml_reporter;

    if (cmd.get_option("--loco")->count() == 0)
    {
        gelex::reml::Estimator estimator(
            cmd.get_option("--max-iter")->as<int>(),
            cmd.get_option("--tol")->as<double>(),
            reml_reporter.as_observer());

        reporter.show_reml_started("");

        auto reml = estimator.fit(model, state);
        reml_reporter.show_result(
            model,
            state,
            estimator.is_converged(),
            estimator.iter_count(),
            cmd.get_option("--max-iter")->as<int>(),
            estimator.loglike());

        auto chr_groups = gelex::build_chr_groups(false, bim);

        reporter.start_scan(
            total_snps, cmd.get_option("--chunk-size")->as<int>(), false);

        scan_groups(chr_groups, reml);
    }
    else
    {
        if (model.random().size() != assoc_data.grm_prefixes.size())
        {
            throw gelex::GelexException(
                "Number of random components in model does not match "
                "number of GRMs provided.");
        }

        std::vector<gelex::LocoReader> loco_readers;
        loco_readers.reserve(assoc_data.grm_prefixes.size());
        for (const auto& path : assoc_data.grm_prefixes)
        {
            loco_readers.emplace_back(path, assoc_data.sample_index);
        }

        auto chr_groups = gelex::build_chr_groups(true, bim);

        reporter.start_scan(
            total_snps, cmd.get_option("--chunk-size")->as<int>(), true);

        std::vector<gelex::LocoRemlResult> loco_results;

        for (const auto& group : chr_groups)
        {
            for (std::size_t i = 0; i < loco_readers.size(); ++i)
            {
                const auto chr_grm_prefix
                    = assoc_data.grm_prefixes[i] + ".chr" + group.name;
                loco_readers[i].load_loco_grm(
                    chr_grm_prefix,
                    assoc_data.sample_index,
                    model.random()[i].K);
            }

            reporter.show_loco_phase(group.name, "REML");

            gelex::reml::Estimator estimator(
                cmd.get_option("--max-iter")->as<int>(),
                cmd.get_option("--tol")->as<double>());
            auto reml = estimator.fit(model, state);

            {
                gelex::LocoRemlResult r;
                r.chr_name = group.name;
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

            reporter.show_loco_phase(group.name, "SCAN");

            scan_groups({group}, reml);
        }

        reporter.show_loco_reml_summary(loco_results);
    }

    reporter.show_complete(cmd.get_option("--out")->as<std::string>());

    return 0;
}
