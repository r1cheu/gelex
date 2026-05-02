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

#include "gelex/model/bayes/checkpoint_reader.h"

#include <fmt/format.h>
#include <cstdint>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "gelex/io/binary_reader.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"

namespace gelex
{

namespace
{

template <typename eT>
auto read_vec(const io::detail::BinaryReader& reader, std::string_view path)
    -> Eigen::VectorX<eT>
{
    auto map = reader.to_map<eT>(path);
    return map.col(0);
}

auto read_scalar(const io::detail::BinaryReader& reader, std::string_view path)
    -> double
{
    auto map = reader.to_map<double>(path);
    return map(0, 0);
}

auto rebuild_pi_count(
    const bayes::TrackerVector& tracker,
    Eigen::Index n_components) -> Eigen::VectorXi
{
    Eigen::VectorXi count = Eigen::VectorXi::Zero(n_components);
    for (Eigen::Index j = 0; j < tracker.size(); ++j)
    {
        const auto comp = tracker(j);
        if (comp >= 0 && comp < n_components)
        {
            ++count(comp);
        }
    }
    return count;
}

auto read_fixed(const io::detail::BinaryReader& reader) -> bayes::FixedState
{
    return bayes::FixedState(
        read_vec<double>(reader, fmt::format("{}/coeff", EffectType::fixed())));
}

auto read_random(const io::detail::BinaryReader& reader)
    -> std::vector<bayes::RandomState>
{
    std::vector<bayes::RandomState> result;
    for (uint8_t i = 0;; ++i)
    {
        const auto coeff_path
            = fmt::format("{}/coeff/{}", EffectType::random(), i);
        if (!reader.contains(coeff_path))
        {
            break;
        }
        result.emplace_back(
            read_vec<double>(reader, coeff_path),
            read_scalar(
                reader,
                fmt::format("{}/variance/{}", EffectType::random(), i)));
    }
    return result;
}

auto read_assignment(
    const io::detail::BinaryReader& reader,
    const EffectType& effect,
    std::string_view group_key) -> bayes::Assignment
{
    auto tracker = read_vec<int8_t>(
        reader, fmt::format("{}/{}/assignment", effect, group_key));
    auto proportion = read_vec<double>(
        reader, fmt::format("{}/{}/proportion", effect, group_key));
    auto count = rebuild_pi_count(tracker, proportion.size());
    return {std::move(tracker), std::move(proportion), std::move(count)};
}

auto read_component_allocation(
    const io::detail::BinaryReader& reader,
    const EffectType& effect) -> bayes::ComponentAllocation
{
    auto assignment = read_assignment(reader, effect, "group");

    std::vector<Eigen::VectorXd> component_u;
    for (uint8_t ci = 0;; ++ci)
    {
        const auto path = fmt::format("{}/component_u/{}", effect, ci);
        if (!reader.contains(path))
        {
            break;
        }
        component_u.push_back(read_vec<double>(reader, path));
    }

    Eigen::VectorXd component_variance;
    if (const auto path = fmt::format("{}/component_variance", effect);
        reader.contains(path))
    {
        component_variance = read_vec<double>(reader, path);
    }

    return {
        std::move(assignment),
        std::move(component_u),
        std::move(component_variance)};
}

auto read_marker_group(
    const io::detail::BinaryReader& reader,
    const EffectType& effect,
    const bayes::MarkerPrior& marker_prior) -> bayes::MarkerAllocation
{
    if (std::holds_alternative<bayes::SpikePrior>(marker_prior))
    {
        return read_assignment(reader, effect, "group");
    }
    return read_component_allocation(reader, effect);
}

auto read_sign(const io::detail::BinaryReader& reader, const EffectType& effect)
    -> bayes::Assignment
{
    return read_assignment(reader, effect, "sign");
}

auto read_one_genetic(
    const io::detail::BinaryReader& reader,
    const bayes::GeneticPrior& prior) -> bayes::GeneticState
{
    const auto effect = EffectType::from_genetic(prior.type);

    auto coeffs = read_vec<double>(reader, fmt::format("{}/coeff", effect));
    auto variance = read_scalar(reader, fmt::format("{}/variance", effect));
    auto u = read_vec<double>(reader, fmt::format("{}/gebv", effect));
    auto marker_variance
        = read_vec<double>(reader, fmt::format("{}/marker_variance", effect));

    std::optional<bayes::MarkerAllocation> group;
    if (!std::holds_alternative<bayes::ContinuousPrior>(prior.marker))
    {
        group = read_marker_group(reader, effect, prior.marker);
    }

    std::optional<bayes::Assignment> sign;
    if (prior.sign)
    {
        sign = read_sign(reader, effect);
    }

    return bayes::GeneticState(
        prior.type,
        std::move(coeffs),
        std::move(u),
        variance,
        std::move(marker_variance),
        std::move(group),
        std::move(sign));
}

auto read_genetics(
    const io::detail::BinaryReader& reader,
    const bayes::Priors& priors) -> std::vector<bayes::GeneticState>
{
    std::vector<bayes::GeneticState> result;
    for (const auto& gp : priors.genetics())
    {
        result.push_back(read_one_genetic(reader, gp));
    }
    return result;
}

auto read_residual(const io::detail::BinaryReader& reader)
    -> bayes::ResidualState
{
    auto y_adj = read_vec<double>(
        reader, fmt::format("{}/adj_pheno", EffectType::residual()));
    auto variance = read_scalar(
        reader, fmt::format("{}/variance", EffectType::residual()));
    return {std::move(y_adj), variance};
}

auto read_rng(const io::detail::BinaryReader& reader) -> std::mt19937_64
{
    auto rng_map = reader.to_map<uint8_t>("rng_state");
    std::string str(
        reinterpret_cast<const char*>(rng_map.data()),
        static_cast<size_t>(rng_map.size()));
    std::istringstream iss(str);
    std::mt19937_64 rng{};
    iss >> rng;
    return rng;
}

auto read_variance_prior(
    const io::detail::BinaryReader& reader,
    std::string_view prefix) -> bayes::VariancePrior
{
    auto nu_map = reader.to_map<double>(fmt::format("{}/nu", prefix));
    auto s2_map = reader.to_map<double>(fmt::format("{}/s2", prefix));
    auto init_map = reader.to_map<double>(fmt::format("{}/init", prefix));
    auto size_map = reader.to_map<double>(fmt::format("{}/size", prefix));
    return {
        .param = {.nu = nu_map(0, 0), .s2 = s2_map(0, 0)},
        .init = init_map(0, 0),
        .size = static_cast<Eigen::Index>(size_map(0, 0))};
}

auto read_genetic_prior(
    const io::detail::BinaryReader& reader,
    const EffectType& effect) -> bayes::GeneticPrior
{
    const auto prefix = fmt::format("prior/{}", effect);
    auto type_map = reader.to_map<uint8_t>(fmt::format("{}/type", prefix));
    const auto type_index = type_map(0, 0);

    auto vp = read_variance_prior(reader, fmt::format("{}/variance", prefix));

    bayes::MarkerPrior marker;
    if (type_index == 0)
    {
        marker = bayes::ContinuousPrior{.variance = vp};
    }
    else if (type_index == 1)
    {
        auto pi_init
            = reader.to_mat<double>(fmt::format("{}/proportion/init", prefix));
        auto est_map = reader.to_map<uint8_t>(
            fmt::format("{}/proportion/estimate", prefix));
        marker = bayes::SpikePrior{
            .variance = vp,
            .proportion
            = {.init = pi_init.col(0), .estimate = est_map(0, 0) != 0}};
    }
    else
    {
        auto pi_init
            = reader.to_mat<double>(fmt::format("{}/proportion/init", prefix));
        auto est_map = reader.to_map<uint8_t>(
            fmt::format("{}/proportion/estimate", prefix));
        auto mult = reader.to_mat<double>(fmt::format("{}/multiplier", prefix));
        marker = bayes::MixturePrior{
            .variance = vp,
            .proportion
            = {.init = pi_init.col(0), .estimate = est_map(0, 0) != 0},
            .multiplier = mult.col(0)};
    }

    std::optional<bayes::SignPrior> sign;
    if (const auto path = fmt::format("{}/sign", prefix); reader.contains(path))
    {
        auto sign_map = reader.to_map<double>(path);
        sign = bayes::SignPrior{.init_value = sign_map(0, 0)};
    }

    return {
        .type = *effect.genetic_mode,
        .marker = std::move(marker),
        .sign = sign};
}

auto read_priors(const io::detail::BinaryReader& reader) -> bayes::Priors
{
    auto residual = read_variance_prior(
        reader, fmt::format("prior/{}", EffectType::residual()));

    std::vector<bayes::RandomPrior> random;
    for (uint8_t i = 0;; ++i)
    {
        const auto path = fmt::format("prior/{}/{}", EffectType::random(), i);
        if (!reader.contains(fmt::format("{}/nu", path)))
        {
            break;
        }
        random.push_back(read_variance_prior(reader, path));
    }

    std::vector<bayes::GeneticPrior> genetics;
    for (auto kind : {GeneticMode::A, GeneticMode::D})
    {
        const auto effect = EffectType::from_genetic(kind);
        const auto path = fmt::format("prior/{}/type", effect);
        if (!reader.contains(path))
        {
            continue;
        }
        genetics.push_back(read_genetic_prior(reader, effect));
    }

    return {std::move(genetics), std::move(random), residual};
}

}  // namespace

auto read_checkpoint(const std::filesystem::path& path) -> Checkpoint
{
    io::detail::BinaryReader reader(path.string());
    auto priors = read_priors(reader);

    auto fixed = read_fixed(reader);
    auto random = read_random(reader);
    auto genetics = read_genetics(reader, priors);
    auto residual = read_residual(reader);
    auto rng = read_rng(reader);

    mcmc::State state(
        std::move(fixed),
        std::move(random),
        std::move(genetics),
        std::move(residual));

    return Checkpoint{std::move(state), rng, std::move(priors)};
}

}  // namespace gelex
