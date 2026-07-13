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

#include <utility>

#include "gelex/algo/reml/estimator.h"
#include "gelex/algo/reml/summary.h"
#include "gelex/exception.h"
#include "gelex/freq/model.h"
#include "gelex/io/reml.h"

#include "cli/cli_helper.h"
#include "cli/common_data.h"
#include "cli/formatter.h"
#include "cli/reml_data.h"
#include "cli/reml_reporter.h"
#include "cli/report_printer.h"

auto reml_execute(const cli::RemlConfig& config) -> int
{
    if (config.random.grm.empty() && !config.random.drand_path
        && config.random.qrand_paths.empty()
        && config.random.interactions.empty())
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

    cli::RemlReporter reml_reporter;
    reml_reporter.show_dataset_summary(model);

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
