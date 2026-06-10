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

#include <iterator>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>

#include "cli/cli_helper.h"
#include "cli/common_data.h"
#include "cli/data_pipe_config.h"
#include "gelex/algo/reml/estimator.h"
#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/freq/design.h"
#include "gelex/freq/model.h"
#include "gelex/io/reml.h"
#include "gelex/types/fixed_designs.h"
#include "reporter.h"

class RemlDataHandler
{
   public:
    auto load_indices(
        CLI::App& cmd,
        std::vector<gelex::dataframe::Index<std::string>*>& indices) -> void
    {
        rand_ = cmd.get_option("--rand")->count() > 0
                    ? std::make_optional(
                          gelex::read_dcovar(
                              cmd.get_option("--rand")->as<std::string>()))
                    : std::nullopt;
        if (rand_)
        {
            indices.push_back(&rand_->index());
        }

        grm_prefixes_ = cmd.get_option("--grm")->as<std::vector<std::string>>();
        grm_indices_.reserve(grm_prefixes_.size());
        for (const auto& path : grm_prefixes_)
        {
            grm_indices_.emplace_back(gelex::read_grm_ids(path));
            indices.push_back(&grm_indices_.back());
        }
    }

    auto gather(const gelex::dataframe::Index<std::string>& common_index)
        -> void
    {
        if (rand_)
        {
            random_designs_ = gelex::make_random_designs(*rand_);
        }
        auto grm_designs = gelex::make_grm_designs(grm_prefixes_, common_index);
        random_designs_.insert(
            random_designs_.end(),
            std::make_move_iterator(grm_designs.begin()),
            std::make_move_iterator(grm_designs.end()));
    }

    auto results() && -> std::vector<gelex::freq::RandomDesign>
    {
        return std::move(random_designs_);
    }

   private:
    std::vector<std::string> grm_prefixes_;
    std::vector<gelex::dataframe::Index<std::string>> grm_indices_;
    std::vector<gelex::freq::RandomDesign> random_designs_;

    std::optional<gelex::dataframe::DataFrame<std::string>> rand_;
};

auto reml_execute(CLI::App& cmd) -> int
{
    int threads = cmd.get_option("--threads")->as<int>();
    cli::setup_parallelization(threads);

    cli::RemlCommandReporter reporter;
    reporter.show_banner();
    reporter.show_config(cmd);

    RemlDataHandler handler;
    cli::BaseData data = cli::load_base_data(handler, cmd);

    gelex::FreqModel model(
        std::move(data.phenotype),
        std::move(data.fixed_design),
        std::move(handler).results());

    gelex::reml::Estimator estimator(
        cmd.get_option("--max-iter")->as<int>(),
        cmd.get_option("--tol")->as<double>(),
        reporter.as_observer());

    gelex::FreqState state(model);
    estimator.fit(model, state);

    const auto out_prefix = cmd.get_option("--out")->as<std::string>();
    gelex::reml::write_summary(model, state, out_prefix);
    gelex::reml::write_effects(model, state, data.sample_ids, out_prefix);

    return 0;
}
