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

#include "reml_data.h"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <iterator>
#include <ranges>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/grm/io.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/freq/design_factory.h"

namespace cli
{

namespace
{

auto split_interaction(const std::string& spec)
    -> std::pair<std::string, std::string>
{
    auto sep = spec.find(':');
    if (sep == std::string::npos || sep == 0 || sep + 1 == spec.size())
    {
        throw gelex::GelexException(
            fmt::format(
                "--interaction expects '<name_a>:<name_b>', got '{}'", spec));
    }
    return {spec.substr(0, sep), spec.substr(sep + 1)};
}

}  // namespace

RemlDataLoader::RemlDataLoader(const RemlDataConfig& config) noexcept
    : config_(config)
{
}

auto RemlDataLoader::load_indices(
    std::vector<const gelex::DataFrameIndex<std::string>*>& indices) -> void
{
    // Names of effects loaded as main components; an interaction operand
    // matching one reuses it instead of reading its GRM from a path.
    std::set<std::string> known_names;

    if (config_.drand_path)
    {
        drand_ = gelex::read_dcovar(*config_.drand_path);
        indices.push_back(&drand_->index());
        for (auto& name : drand_->names())
        {
            known_names.insert(std::move(name));
        }
    }

    qrand_.reserve(config_.qrand_paths.size());
    for (const auto& path : config_.qrand_paths)
    {
        qrand_.emplace_back(gelex::read_qcovar(path));
        indices.push_back(&qrand_.back().index());
        known_names.insert(std::filesystem::path(path).stem().string());
    }

    grm_indices_.reserve(config_.grm.size());
    for (const auto& path : config_.grm)
    {
        grm_indices_.emplace_back(gelex::read_grm_ids(path));
        indices.push_back(&grm_indices_.back());
        known_names.insert(std::filesystem::path(path).filename().string());
    }

    // Load an interaction operand's GRM on first sight. Operands naming an
    // already-loaded effect or an already-registered path GRM are skipped (both
    // are resolved by name in gather); any other operand is treated as a GRM
    // prefix whose ids are read here so it joins the sample intersection, with
    // its matrix read later in gather.
    auto load_operand_grm = [&](const std::string& operand)
    {
        if (known_names.contains(operand)
            || interaction_grms_.contains(operand))
        {
            return;
        }
        auto [it, _] = interaction_grms_.emplace(
            operand,
            InteractionGrm{.index = gelex::read_grm_ids(operand), .K = {}});
        indices.push_back(&it->second.index);
    };

    for (const auto& spec : config_.interactions)
    {
        auto [lhs, rhs] = split_interaction(spec);
        load_operand_grm(lhs);
        load_operand_grm(rhs);
    }
}

auto RemlDataLoader::gather(
    const gelex::DataFrameIndex<std::string>& common_index) -> void
{
    if (drand_)
    {
        drand_->gather(common_index);
        random_designs_ = gelex::freq::make_random_designs(*drand_);
    }
    for (auto&& [i, frame] : std::views::enumerate(qrand_))
    {
        frame.gather(common_index);
        random_designs_.push_back(
            gelex::freq::make_quantitative_random_design(
                frame,
                std::filesystem::path(
                    config_.qrand_paths[static_cast<std::size_t>(i)])
                    .stem()
                    .string()));
    }
    auto grm_designs = gelex::freq::make_grm_designs(config_.grm, common_index);
    random_designs_.insert(
        random_designs_.end(),
        std::make_move_iterator(grm_designs.begin()),
        std::make_move_iterator(grm_designs.end()));

    for (auto& [prefix, grm] : interaction_grms_)
    {
        grm.K = gelex::read_grm(prefix, &common_index);
    }

    // Both kernels are read fully before push_back, so the references
    // resolve_operand returns never outlive the make_interaction_design call.
    for (const auto& spec : config_.interactions)
    {
        auto [lhs, rhs] = split_interaction(spec);
        random_designs_.push_back(
            gelex::freq::make_interaction_design(
                spec, resolve_operand(lhs), resolve_operand(rhs)));
    }
}

auto RemlDataLoader::resolve_operand(const std::string& name) const
    -> const Eigen::MatrixXd&
{
    if (auto it = interaction_grms_.find(name); it != interaction_grms_.end())
    {
        return it->second.K;
    }
    auto it = std::ranges::find(
        random_designs_, name, &gelex::freq::RandomDesign::name);
    if (it == random_designs_.end())
    {
        auto names = random_designs_
                     | std::views::transform(&gelex::freq::RandomDesign::name);
        throw gelex::GelexException(
            fmt::format(
                "--interaction references unknown effect '{}'; loaded effects "
                "are [{}]",
                name,
                fmt::join(names, ", ")));
    }
    return it->K;
}

auto RemlDataLoader::results() && -> std::vector<gelex::freq::RandomDesign>
{
    return std::move(random_designs_);
}

}  // namespace cli
