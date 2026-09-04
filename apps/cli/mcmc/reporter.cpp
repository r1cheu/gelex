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
#include <fmt/color.h>
#include <fmt/format.h>
#include <string>
#include <string_view>

#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/progress.h"
#include "gelex/bayes/mcmc_progress.h"
#include "gelex/bayes/model.h"
#include "gelex/genetic_mode.h"

#include "cli/formatter.h"
#include "cli/report_printer.h"
#include "cli/theme.h"

namespace cli
{
namespace
{

auto phase_prefix(std::string_view phase, ColorRole color) -> std::string
{
    return fmt::format(
        "[{}]",
        colorize(style_for(color) | fmt::emphasis::bold, std::string{phase}));
}

}  // namespace

GenoReporter::GenoReporter(std::size_t total)
    : progress_{"", total, "SNP"},
      estimate_rate_{cli::make_rate()},
      estimate_eta_{cli::make_eta(total)}
{
}

auto GenoReporter::show_total(int64_t num_variants) -> void
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
        progress_.finish();
        cli::printer().on_progress_finished();
        return;
    }

    progress_.update(
        {.current = event.current,
         .rate = estimate_rate_(event.current),
         .eta = estimate_eta_(event.current)});
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

auto McmcReporter::show_sampling_header() -> void
{
    cli::printer().block(cli::section("MCMC Sampling:"));
}

McmcReporter::McmcReporter(std::size_t total, std::size_t burn_in)
    : burn_in_{burn_in},
      progress_{
          burn_in == 0 ? phase_prefix("SAMPLE", ColorRole::success)
                       : phase_prefix("BURN-IN", ColorRole::warning),
          total,
          "iter"},
      estimate_rate_{cli::make_rate()},
      estimate_eta_{cli::make_eta(total)}
{
}

auto McmcReporter::on_event(const gelex::MCMCProgressEvent& event) -> void
{
    if (event.done)
    {
        progress_.finish();
        cli::printer().on_progress_finished();
        return;
    }

    progress_.update(
        {.current = event.current,
         .rate = estimate_rate_(event.current),
         .eta = estimate_eta_(event.current)});

    if (event.current == burn_in_)
    {
        progress_.set_prefix(phase_prefix("SAMPLE", ColorRole::success));
    }
}

}  // namespace cli
