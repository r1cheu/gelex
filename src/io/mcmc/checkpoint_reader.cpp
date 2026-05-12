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

#include "gelex/io/mcmc/checkpoint_reader.h"

#include <fmt/format.h>
#include <cstdint>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/io/detail/binary_reader.h"
#include "gelex/model/bayes/legacy_method.h"
#include "gelex/model/bayes/prior.h"

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

auto read_uint8(const io::detail::BinaryReader& reader, std::string_view path)
    -> uint8_t
{
    auto map = reader.to_map<uint8_t>(path);
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

auto read_sign(const io::detail::BinaryReader& reader, const EffectType& effect)
    -> bayes::Assignment
{
    return read_assignment(reader, effect, "sign");
}

auto read_one_genetic(
    const io::detail::BinaryReader& reader,
    GeneticMode mode,
    bool has_mixture,
    bool has_sign) -> bayes::GeneticState
{
    const auto effect = EffectType::from_genetic(mode);

    auto coeffs = read_vec<double>(reader, fmt::format("{}/coeff", effect));
    auto variance = read_scalar(reader, fmt::format("{}/variance", effect));
    auto u = read_vec<double>(reader, fmt::format("{}/gebv", effect));
    auto marker_variance
        = read_vec<double>(reader, fmt::format("{}/marker_variance", effect));

    std::optional<bayes::MarkerAllocation> group;
    if (has_mixture)
    {
        if (reader.contains(fmt::format("{}/component_u/0", effect)))
        {
            group = read_component_allocation(reader, effect);
        }
        else
        {
            group = read_assignment(reader, effect, "group");
        }
    }

    std::optional<bayes::Assignment> sign;
    if (has_sign)
    {
        sign = read_sign(reader, effect);
    }

    return bayes::GeneticState(
        mode,
        std::move(coeffs),
        std::move(u),
        variance,
        std::move(marker_variance),
        std::move(group),
        std::move(sign));
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

auto read_variance_spec(
    const io::detail::BinaryReader& reader,
    std::string_view prefix) -> bayes::OldVarianceSpec
{
    const auto nu = read_scalar(reader, fmt::format("{}/nu", prefix));
    const auto s2 = read_scalar(reader, fmt::format("{}/s2", prefix));
    const auto init = read_scalar(reader, fmt::format("{}/init", prefix));
    const auto scope_raw = read_uint8(reader, fmt::format("{}/scope", prefix));
    return {
        .scope = static_cast<bayes::MarkerVarianceScope>(scope_raw),
        .init = init,
        .prior = {nu, s2},
    };
}

auto read_categorical_spec(
    const io::detail::BinaryReader& reader,
    std::string_view prefix) -> bayes::CategoricalSpec
{
    auto init = read_vec<double>(reader, fmt::format("{}/init", prefix));
    const auto estimate
        = read_uint8(reader, fmt::format("{}/estimate", prefix));
    const Eigen::Index n = init.size();
    return {
        .init = std::move(init),
        .prior = {Eigen::VectorXi::Ones(n)},
        .estimate = estimate != 0,
    };
}

auto read_genetic_spec(
    const io::detail::BinaryReader& reader,
    std::string_view prefix) -> bayes::GeneticSpec
{
    const auto mode_raw = read_uint8(reader, fmt::format("{}/mode", prefix));
    auto variance
        = read_variance_spec(reader, fmt::format("{}/variance", prefix));

    std::optional<bayes::CategoricalSpec> sign;
    if (reader.contains(fmt::format("{}/sign/init", prefix)))
    {
        sign = read_categorical_spec(reader, fmt::format("{}/sign", prefix));
    }

    return {
        .mode = static_cast<GeneticMode>(mode_raw),
        .variance = std::move(variance),
        .sign = std::move(sign),
    };
}

auto read_genetic_prior(
    const io::detail::BinaryReader& reader,
    std::string_view prefix) -> bayes::OldGeneticPrior
{
    const auto kind = read_uint8(reader, fmt::format("{}/kind", prefix));

    std::variant<bayes::GeneticSpec, bayes::JointSpec> spec;
    if (kind == 0)
    {
        spec = read_genetic_spec(reader, fmt::format("{}/spec", prefix));
    }
    else
    {
        spec = bayes::JointSpec{
            read_genetic_spec(reader, fmt::format("{}/spec/additive", prefix)),
            read_genetic_spec(reader, fmt::format("{}/spec/dominance", prefix)),
        };
    }

    std::optional<bayes::Mixture> mixture;
    if (read_uint8(reader, fmt::format("{}/mixture/present", prefix)) != 0)
    {
        const auto strategy_idx
            = read_uint8(reader, fmt::format("{}/mixture/strategy", prefix));
        auto proportions = read_categorical_spec(
            reader, fmt::format("{}/mixture/proportions", prefix));

        bayes::VarianceStrategy strategy;
        if (strategy_idx == 0)
        {
            strategy = bayes::SpikeSlab{};
        }
        else if (strategy_idx == 1)
        {
            auto multiplier = read_vec<double>(
                reader, fmt::format("{}/mixture/multiplier", prefix));
            strategy = bayes::ScaledMixture{std::move(multiplier)};
        }
        else
        {
            strategy = bayes::JointMixture{};
        }

        mixture = bayes::Mixture{std::move(strategy), std::move(proportions)};
    }

    return {std::move(spec), std::move(mixture)};
}

auto read_method(const io::detail::BinaryReader& reader)
    -> bayes::LegacyBayesMethod
{
    bayes::LegacyBayesMethod method;

    for (uint8_t i = 0;; ++i)
    {
        const auto prefix = fmt::format("method/genetic/{}", i);
        if (!reader.contains(fmt::format("{}/kind", prefix)))
        {
            break;
        }
        method.genetics.push_back(read_genetic_prior(reader, prefix));
    }

    for (uint8_t i = 0;; ++i)
    {
        const auto prefix = fmt::format("method/random/{}", i);
        if (!reader.contains(fmt::format("{}/nu", prefix)))
        {
            break;
        }
        method.randoms.push_back(read_variance_spec(reader, prefix));
    }

    method.residual = read_variance_spec(reader, "method/residual");
    return method;
}

struct GeneticModeInfo
{
    GeneticMode mode;
    bool has_mixture;
    bool has_sign;
};

auto collect_genetic_modes(const bayes::LegacyBayesMethod& method)
    -> std::vector<GeneticModeInfo>
{
    std::vector<GeneticModeInfo> result;
    for (const auto& prior : method.genetics)
    {
        const bool has_mixture = prior.mixture.has_value();
        std::visit(
            [&](const auto& spec)
            {
                using T = std::decay_t<decltype(spec)>;
                if constexpr (std::is_same_v<T, bayes::GeneticSpec>)
                {
                    result.push_back(
                        {spec.mode, has_mixture, spec.sign.has_value()});
                }
                else
                {
                    result.push_back(
                        {spec.additive.mode,
                         has_mixture,
                         spec.additive.sign.has_value()});
                    result.push_back(
                        {spec.dominance.mode,
                         has_mixture,
                         spec.dominance.sign.has_value()});
                }
            },
            prior.spec);
    }
    return result;
}

auto read_genetics(
    const io::detail::BinaryReader& reader,
    const bayes::LegacyBayesMethod& method) -> std::vector<bayes::GeneticState>
{
    const auto modes = collect_genetic_modes(method);
    std::vector<bayes::GeneticState> result;
    result.reserve(modes.size());
    for (const auto& [mode, has_mixture, has_sign] : modes)  // NOLINT
    {
        result.push_back(read_one_genetic(reader, mode, has_mixture, has_sign));
    }
    return result;
}

}  // namespace

auto read_checkpoint(const std::filesystem::path& path) -> Checkpoint
{
    io::detail::BinaryReader reader(path.string());

    auto method = read_method(reader);

    auto fixed = read_fixed(reader);
    auto random = read_random(reader);
    auto genetics = read_genetics(reader, method);
    auto residual = read_residual(reader);
    auto rng = read_rng(reader);

    mcmc::State state(
        std::move(fixed),
        std::move(random),
        std::move(genetics),
        std::move(residual));

    return Checkpoint{std::move(state), rng, std::move(method)};
}

}  // namespace gelex
