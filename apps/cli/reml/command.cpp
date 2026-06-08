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

#include <argparse.h>
#include <optional>

#include "cli/cli_helper.h"
#include "cli/data_pipe_config.h"
#include "cli/reml_reporter.h"
#include "gelex/algo/reml/estimator.h"
#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/freq/design.h"
#include "gelex/freq/model.h"
#include "gelex/types/fixed_designs.h"

auto reml_execute(argparse::ArgumentParser& cmd) -> int
{
    int threads = cmd.get<int>("--threads");
    gelex::cli::setup_parallelization(threads);

    gelex::cli::RemlReporter reporter;

    std::vector<gelex::dataframe::Index<std::string>*> indices;

    // read dataset
    auto pheno_col = cmd.get<int>("--pheno-col");
    if (pheno_col < 2)
    {
        throw gelex::GelexException("--pheno-col must be >= 2");
    }
    auto pheno_col_offset = static_cast<std::size_t>(pheno_col - 2);
    auto phenotype = gelex::read_pheno(cmd.get("--pheno"), &pheno_col_offset);
    std::cout << phenotype.col(0).name() << " is used as the phenotype.\n";
    indices.push_back(&phenotype.index());

    std::optional<gelex::dataframe::DataFrame<std::string>> qcovar;
    std::optional<gelex::dataframe::DataFrame<std::string>> dcovar;
    std::optional<gelex::dataframe::DataFrame<std::string>> rand;
    if (cmd.is_used("--qcovar"))
    {
        qcovar = std::make_optional(gelex::read_qcovar(cmd.get("--qcovar")));
        indices.push_back(&qcovar->index());
    }
    if (cmd.is_used("--dcovar"))
    {
        dcovar = std::make_optional(gelex::read_dcovar(cmd.get("--dcovar")));
        indices.push_back(&dcovar->index());
    }
    if (cmd.is_used("--rand"))
    {
        rand = std::make_optional(gelex::read_dcovar(cmd.get("--rand")));
        indices.push_back(&rand->index());
    }

    auto grm_prefixes = cmd.get<std::vector<std::string>>("--grm");
    std::vector<gelex::dataframe::Index<std::string>> grm_indices;
    grm_indices.reserve(grm_prefixes.size());
    for (const auto& path : grm_prefixes)
    {
        grm_indices.emplace_back(gelex::read_grm_ids(path));
        indices.push_back(&grm_indices.back());
    }

    auto common_index = gelex::dataframe::intersect<std::string>(indices);

    phenotype.gather(common_index);
    if (qcovar)
    {
        qcovar->gather(common_index);
    }
    if (dcovar)
    {
        dcovar->gather(common_index);
    }
    if (rand)
    {
        rand->gather(common_index);
    }
    for (auto& idx : grm_indices)
    {
        idx.gather(common_index);
    }

    std::optional<gelex::FixedDesign> fixed_design;

    if (!qcovar && !dcovar)
    {
        fixed_design
            = std::make_optional(gelex::FixedDesign::make(common_index.size()));
    }
    else
    {
        fixed_design = std::make_optional(
            gelex::FixedDesign::make(
                qcovar ? std::make_optional(
                             gelex::make_quantitative_covariate(*qcovar))
                       : std::nullopt,
                dcovar ? std::make_optional(
                             gelex::make_discrete_covariate(*dcovar))
                       : std::nullopt));
    }

    auto random_designs = rand ? gelex::make_random_designs(*rand)
                               : std::vector<gelex::freq::RandomDesign>{};
    auto genetic_designs
        = gelex::make_genetic_designs(grm_prefixes, common_index);

    random_designs.insert(
        random_designs.end(),
        std::make_move_iterator(genetic_designs.begin()),
        std::make_move_iterator(genetic_designs.end()));

    gelex::FreqModel model(
        phenotype.col(0).to_mat<double>(),
        std::move(*fixed_design),
        std::move(random_designs),
        {});

    gelex::reml::Estimator estimator(
        cmd.get<int>("--max-iter"),
        cmd.get<double>("--tol"),
        reporter.as_observer());

    gelex::FreqState state(model);
    estimator.fit(model, state);

    return 0;
}
