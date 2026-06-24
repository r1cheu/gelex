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

#include "cli/cli_helper.h"
#include "cli/common_data.h"
#include "cli/reml_reporter.h"
#include "gelex/algo/reml/estimator.h"
#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/freq/design.h"
#include "gelex/freq/model.h"
#include "gelex/io/reml.h"
#include "gelex/types/fixed_designs.h"

class RemlDataHandler
{
   public:
    explicit RemlDataHandler(const cli::RemlConfig& config) noexcept
        : config_(config)
    {
    }

    auto load_indices(
        std::vector<gelex::dataframe::Index<std::string>*>& indices) -> void
    {
        rand_ = config_.rand_path
                    ? std::make_optional(gelex::read_dcovar(*config_.rand_path))
                    : std::nullopt;
        if (rand_)
        {
            indices.push_back(&rand_->index());
        }

        grm_indices_.reserve(config_.grm_prefixes.size());
        for (const auto& path : config_.grm_prefixes)
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
            rand_->gather(common_index);
            random_designs_ = gelex::make_random_designs(*rand_);
        }
        auto grm_designs
            = gelex::make_grm_designs(config_.grm_prefixes, common_index);
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
    const cli::RemlConfig& config_;
    std::vector<gelex::dataframe::Index<std::string>> grm_indices_;
    std::vector<gelex::freq::RandomDesign> random_designs_;

    std::optional<gelex::dataframe::DataFrame<std::string>> rand_;
};

auto reml_execute(const cli::RemlConfig& config) -> int
{
    cli::setup_parallelization(config.threads);

    RemlDataHandler handler(config);
    cli::BaseData data = cli::load_base_data(handler, config.base_data);

    gelex::FreqModel model(
        std::move(data.phenotype),
        std::move(data.fixed_design),
        std::move(handler).results());

    cli::RemlReporter reml_reporter;
    gelex::reml::Estimator estimator(
        config.max_iter, config.tolerance, reml_reporter.as_observer());

    gelex::FreqState state(model);
    estimator.fit(model, state);
    reml_reporter.show_result(
        model,
        state,
        estimator.is_converged(),
        estimator.iter_count(),
        config.max_iter,
        estimator.loglike());

    gelex::reml::write_summary(model, state, config.out_prefix);
    gelex::reml::write_effects(
        model, state, data.sample_ids, config.out_prefix);

    return 0;
}
