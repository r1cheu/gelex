// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "command.h"

#include <fmt/ranges.h>
#include <ranges>
#include <utility>

#include "gelex/exception.h"
#include "gelex/freq/model.h"
#include "gelex/freq/reml/estimator.h"
#include "gelex/freq/reml/io.h"
#include "gelex/freq/reml/summary.h"

#include "cli/common_data.h"
#include "cli/formatter.h"
#include "cli/reml_data.h"
#include "cli/reml_reporter.h"
#include "cli/report_printer.h"
#include "cli/runtime.h"
#include "cli/summary.h"

auto reml_execute(const cli::RemlConfig& config) -> int
{
    if (!config.random.has_random_effect())
    {
        throw gelex::GelexException(
            "REML needs at least one random effect; provide --grm, --drand, "
            "--qrand, or --interaction");
    }

    cli::setup_parallelization(config.threads);

    cli::RemlDataLoader loader(config.random);
    cli::BaseData data = cli::load_base_data(loader, config.base_data);
    auto random_designs = std::move(loader).results();

    gelex::FreqModel model(
        std::move(data.phenotype),
        std::move(data.fixed_design),
        std::move(random_designs));

    cli::Summary{"Dataset Summary"}
        .field("Trait", "{}", data.pheno_name)
        .field("Samples", "{}", model.num_individuals())
        .show();

    const auto random_effect_names
        = model.random()
          | std::views::transform([](const auto& design)
                                  { return design.name; });
    cli::Summary{"Model Summary"}
        .field("Fixed terms", "{}", model.fixed().column_names().size())
        .field("Random effects", "{}", fmt::join(random_effect_names, ", "))
        .show();

    cli::RemlReporter reml_reporter;

    gelex::FreqState state(model);
    gelex::Estimator estimator(
        config.max_iter, config.tolerance, reml_reporter.as_observer());

    auto fit = estimator.fit(model, state);
    reml_reporter.show_result(model, fit.summary, config.max_iter);

    gelex::write_summary(model, state, fit.summary.loglike, config.out_prefix);
    gelex::write_effects(model, state, data.sample_ids, config.out_prefix);
    cli::printer().block(
        cli::results_saved(config.out_prefix, ".summary, .effects, .log"));

    return 0;
}
