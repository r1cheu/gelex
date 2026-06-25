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
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "cli/cli_helper.h"
#include "cli/common_data.h"
#include "cli/reml_reporter.h"
#include "cli/report_printer.h"
#include "gelex/algo/reml/estimator.h"
#include "gelex/algo/reml/loco_result.h"
#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/freq/design.h"
#include "gelex/freq/model.h"
#include "gelex/io/grm/loco_reader.h"
#include "gelex/io/reml.h"
#include "gelex/types/fixed_designs.h"

namespace
{

auto discover_loco_chromosomes(const std::string& grm_prefix)
    -> std::vector<std::string>
{
    std::filesystem::path prefix_path(grm_prefix);
    std::filesystem::path dir = prefix_path.has_parent_path()
                                    ? prefix_path.parent_path()
                                    : std::filesystem::path(".");
    const std::string marker = prefix_path.filename().string() + ".chr";

    std::vector<std::string> chromosomes;
    if (!std::filesystem::exists(dir))
    {
        return chromosomes;
    }
    for (const auto& entry : std::filesystem::directory_iterator(dir))
    {
        const std::string name = entry.path().filename().string();
        if (name.starts_with(marker) && name.ends_with(".bin"))
        {
            chromosomes.push_back(
                name.substr(marker.size(), name.size() - marker.size() - 4));
        }
    }

    const auto as_number
        = [](const std::string& s) -> std::optional<std::int64_t>
    {
        std::int64_t value = 0;
        const auto* end = s.data() + s.size();
        auto [ptr, ec] = std::from_chars(s.data(), end, value);
        if (ec == std::errc{} && ptr == end)
        {
            return value;
        }
        return std::nullopt;
    };
    std::sort(
        chromosomes.begin(),
        chromosomes.end(),
        [&](const std::string& a, const std::string& b)
        {
            const auto na = as_number(a);
            const auto nb = as_number(b);
            if (na && nb)
            {
                return *na < *nb;
            }
            if (na != nb)
            {
                return na.has_value();
            }
            return a < b;
        });
    return chromosomes;
}

}  // namespace

struct RemlData
{
    gelex::dataframe::Index<std::string> sample_index;
    std::vector<gelex::freq::RandomDesign> random_designs;
};

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
        sample_index_ = common_index;
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

    auto results() && -> RemlData
    {
        return RemlData{
            .sample_index = std::move(sample_index_),
            .random_designs = std::move(random_designs_)};
    }

   private:
    const cli::RemlConfig& config_;
    std::vector<gelex::dataframe::Index<std::string>> grm_indices_;
    gelex::dataframe::Index<std::string> sample_index_;
    std::vector<gelex::freq::RandomDesign> random_designs_;

    std::optional<gelex::dataframe::DataFrame<std::string>> rand_;
};

auto reml_execute(const cli::RemlConfig& config) -> int
{
    cli::setup_parallelization(config.threads);

    RemlDataHandler handler(config);
    cli::BaseData data = cli::load_base_data(handler, config.base_data);
    auto reml_data = std::move(handler).results();

    gelex::FreqModel model(
        std::move(data.phenotype),
        std::move(data.fixed_design),
        std::move(reml_data.random_designs));

    gelex::FreqState state(model);

    if (!config.loco)
    {
        cli::RemlReporter reml_reporter;
        gelex::reml::Estimator estimator(
            config.max_iter, config.tolerance, reml_reporter.as_observer());

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

    if (model.random().size() < config.grm_prefixes.size())
    {
        throw gelex::GelexException(
            "Number of random components in model is smaller than the number "
            "of GRMs provided.");
    }
    const std::size_t grm_offset
        = model.random().size() - config.grm_prefixes.size();

    std::vector<gelex::LocoReader> loco_readers;
    loco_readers.reserve(config.grm_prefixes.size());
    for (const auto& path : config.grm_prefixes)
    {
        loco_readers.emplace_back(path, reml_data.sample_index);
    }

    auto chr_names = discover_loco_chromosomes(config.grm_prefixes.front());
    if (chr_names.empty())
    {
        throw gelex::GelexException(
            fmt::format(
                "LOCO error: no chromosome GRM files matching '{}.chr*.bin' "
                "found.",
                config.grm_prefixes.front()));
    }

    cli::printer().line(
        "   Running LOCO REML over {} chromosomes", chr_names.size());

    std::vector<gelex::LocoRemlResult> loco_results;
    loco_results.reserve(chr_names.size());

    for (const auto& chr : chr_names)
    {
        for (std::size_t i = 0; i < loco_readers.size(); ++i)
        {
            const auto chr_grm_prefix = config.grm_prefixes[i] + ".chr" + chr;
            loco_readers[i].load_loco_grm(
                chr_grm_prefix,
                reml_data.sample_index,
                model.random()[grm_offset + i].K);
        }

        gelex::reml::Estimator estimator(config.max_iter, config.tolerance);
        estimator.fit(model, state);

        gelex::LocoRemlResult result;
        result.chr_name = chr;
        result.loglike = estimator.loglike();
        result.converged = estimator.is_converged();
        result.residual_variance = state.residual().variance;
        result.residual_variance_se = state.residual().variance_se;
        for (std::size_t i = 0; i < state.random().size(); ++i)
        {
            const auto& random = state.random()[i];
            result.random.push_back(
                {.name = model.random()[i].name,
                 .variance = random.variance,
                 .variance_se = random.variance_se,
                 .variance_ratio = random.variance_ratio,
                 .variance_ratio_se = random.variance_ratio_se});
        }
        loco_results.push_back(std::move(result));
    }

    cli::print_loco_reml_summary(loco_results);
    gelex::reml::write_loco_summary(loco_results, config.out_prefix);

    return 0;
}
