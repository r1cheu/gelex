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

#include "reporter.h"

#include <Eigen/Core>
#include <cstdint>
#include <fmt/format.h>
#include <string>
#include <string_view>

#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/progress.h"
#include "gelex/bayes/mcmc_progress.h"
#include "gelex/bayes/model.h"
#include "gelex/genetic_mode.h"

#include "cli/formatter.h"
#include "cli/progress_bar.h"
#include "cli/report_printer.h"

namespace cli
{

GenoReporter::GenoReporter() : eta_(1) {}

auto GenoReporter::show_total(int64_t num_variants) const -> void
{
    cli::printer().line(
        "{}", cli::field("Variants", "{}", cli::AbbrNumber(num_variants)));
}

auto GenoReporter::show_loaded(const gelex::bayes::GeneticDesign& design) const
    -> void
{
    for (const gelex::GeneticMode mode : design.each_mode())
    {
        show_loaded_mode(
            mode,
            design.cols(),
            design.cols()
                - static_cast<Eigen::Index>(
                    design.projection(mode).valid_indices().size()));
    }
}

auto GenoReporter::show_loaded_mode(
    gelex::GeneticMode mode,
    int64_t num_snps,
    int64_t invalid_snps) const -> void
{
    const std::string label
        = (mode == gelex::GeneticMode::D) ? "Dominance" : "Additive";
    if (invalid_snps == 0)
    {
        cli::printer().line("{}", cli::field(label, "all valid"));
        return;
    }
    cli::printer().line(
        "{}",
        cli::field(
            label,
            "{} valid ({} excluded)",
            cli::AbbrNumber(num_snps - invalid_snps),
            cli::AbbrNumber(invalid_snps)));
}

auto GenoReporter::on_event(const gelex::GenotypeProgressEvent& event) -> void
{
    if (event.done)
    {
        if (bar_active_)
        {
            bar_.display->done();
            bar_active_ = false;
            cli::clear_finished_line();
            cli::printer().on_progress_finished();
        }
        return;
    }

    if (!bar_active_)
    {
        total_ = event.total;
        progress_ = 0;
        eta_.reset(total_);
        bar_ = cli::create_progress_bar(progress_, total_);
        bar_.display->show();
        bar_active_ = true;
    }

    progress_ = event.current;
    if (bar_.after_bar)
    {
        bar_.after_bar->message(
            fmt::format(
                "{:.1f}% ({}/{} SNPs) | ETA: {}",
                static_cast<double>(progress_) / static_cast<double>(total_)
                    * 100.0,
                cli::AbbrNumber(progress_),
                cli::AbbrNumber(total_),
                eta_.get_eta(progress_)));
    }
}

auto McmcReporter::show_dataset_summary(
    const gelex::BayesModel& model,
    std::string_view pheno_name) -> void
{
    auto& p = cli::printer();
    p.block(cli::section("Dataset Summary:"));
    p.line(cli::field("Trait", "{}", pheno_name));
    p.line(cli::field("Analyzed Samples", "{}", model.num_individuals()));
    p.line(cli::field("Covariates", "{}", model.fixed().X().cols()));
}

auto McmcReporter::on_event(const gelex::MCMCProgressEvent& event) -> void
{
    if (!init_progress_)
    {
        init_progress_ = true;
        cli::printer().block(cli::section("MCMC Sampling:"));
        bar_ = cli::create_progress_bar(
            iter_, event.total, "{bar} {value}/{total} [{speed:.1f}/s]");
        bar_.display->show();
    }

    if (event.done)
    {
        bar_.display->done();
        cli::printer().on_progress_finished();
        return;
    }

    iter_ = event.current;
    bar_.before_bar->message("");
}

}  // namespace cli
