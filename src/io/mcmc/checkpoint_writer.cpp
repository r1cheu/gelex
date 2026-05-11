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

#include "gelex/io/mcmc/checkpoint_writer.h"

#include <fmt/format.h>
#include <cstdint>
#include <span>
#include <sstream>
#include <variant>

#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/io/detail/binary_writer.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"

namespace gelex
{

namespace
{

auto write_scalar(
    io::detail::BinaryWriter& writer,
    std::string_view path,
    double value) -> void
{
    auto handle = writer.reserve<double>(path, 1, 1);
    writer.write(handle, value);
}

auto write_uint8(
    io::detail::BinaryWriter& writer,
    std::string_view path,
    uint8_t value) -> void
{
    auto handle = writer.reserve<uint8_t>(path, 1, 1);
    writer.write(handle, value);
}

auto write_fixed(io::detail::BinaryWriter& writer, const bayes::FixedState& fs)
    -> void
{
    writer.write(fmt::format("{}/coeff", EffectType::fixed()), fs.coeffs);
}

auto write_random(
    io::detail::BinaryWriter& writer,
    const std::vector<bayes::RandomState>& random_states) -> void
{
    for (uint8_t i = 0; i < static_cast<uint8_t>(random_states.size()); ++i)
    {
        const auto& rs = random_states[i];
        writer.write(
            fmt::format("{}/coeff/{}", EffectType::random(), i), rs.coeffs);
        write_scalar(
            writer,
            fmt::format("{}/variance/{}", EffectType::random(), i),
            rs.variance);
    }
}

auto write_assignment(
    io::detail::BinaryWriter& writer,
    const EffectType& effect,
    const bayes::Assignment& a) -> void
{
    writer.write(fmt::format("{}/group/assignment", effect), a.tracker);
    writer.write(fmt::format("{}/group/proportion", effect), a.proportion);
}

auto write_component_allocation(
    io::detail::BinaryWriter& writer,
    const EffectType& effect,
    const bayes::ComponentAllocation& ca) -> void
{
    write_assignment(writer, effect, ca.assignment);

    for (uint8_t ci = 0; ci < static_cast<uint8_t>(ca.component_u.size()); ++ci)
    {
        writer.write(
            fmt::format("{}/component_u/{}", effect, ci), ca.component_u[ci]);
    }

    if (ca.component_variance.size() > 0)
    {
        writer.write(
            fmt::format("{}/component_variance", effect),
            ca.component_variance);
    }
}

auto write_marker_group(
    io::detail::BinaryWriter& writer,
    const EffectType& effect,
    const bayes::MarkerAllocation& group) -> void
{
    std::visit(
        [&](const auto& alloc)
        {
            using T = std::decay_t<decltype(alloc)>;
            if constexpr (std::is_same_v<T, bayes::Assignment>)
            {
                write_assignment(writer, effect, alloc);
            }
            else if constexpr (std::is_same_v<T, bayes::ComponentAllocation>)
            {
                write_component_allocation(writer, effect, alloc);
            }
        },
        group);
}

auto write_genetic_sign(
    io::detail::BinaryWriter& writer,
    const EffectType& effect,
    const bayes::Assignment& sign) -> void
{
    writer.write(fmt::format("{}/sign/assignment", effect), sign.tracker);
    writer.write(fmt::format("{}/sign/proportion", effect), sign.proportion);
}

auto write_genetics(
    io::detail::BinaryWriter& writer,
    const std::vector<bayes::GeneticState>& genetic_states) -> void
{
    for (const auto& gs : genetic_states)
    {
        const auto effect = EffectType::from_genetic(gs.type);
        writer.write(fmt::format("{}/coeff", effect), gs.coeffs);
        write_scalar(writer, fmt::format("{}/variance", effect), gs.variance);
        writer.write(fmt::format("{}/gebv", effect), gs.u);
        writer.write(
            fmt::format("{}/marker_variance", effect), gs.marker_variance);

        if (gs.group)
        {
            write_marker_group(writer, effect, *gs.group);
        }
        if (gs.sign)
        {
            write_genetic_sign(writer, effect, *gs.sign);
        }
    }
}

auto write_residual(
    io::detail::BinaryWriter& writer,
    const bayes::ResidualState& rs) -> void
{
    writer.write(fmt::format("{}/adj_pheno", EffectType::residual()), rs.y_adj);
    write_scalar(
        writer,
        fmt::format("{}/variance", EffectType::residual()),
        rs.variance);
}

auto write_rng(io::detail::BinaryWriter& writer, const std::mt19937_64& rng)
    -> void
{
    std::ostringstream oss;
    oss << rng;
    const auto str = oss.str();
    auto handle = writer.reserve<uint8_t>(
        "rng_state", static_cast<Eigen::Index>(str.size()), 1);
    writer.write(
        handle,
        std::span<const uint8_t>(
            reinterpret_cast<const uint8_t*>(str.data()), str.size()));
}

auto write_variance_spec(
    io::detail::BinaryWriter& writer,
    std::string_view prefix,
    const bayes::OldVarianceSpec& spec) -> void
{
    write_scalar(
        writer, fmt::format("{}/nu", prefix), spec.prior.degrees_of_freedom);
    write_scalar(writer, fmt::format("{}/s2", prefix), spec.prior.scale);
    write_scalar(writer, fmt::format("{}/init", prefix), spec.init);
    write_uint8(
        writer,
        fmt::format("{}/scope", prefix),
        static_cast<uint8_t>(spec.scope));
}

auto write_categorical_spec(
    io::detail::BinaryWriter& writer,
    std::string_view prefix,
    const bayes::CategoricalSpec& spec) -> void
{
    writer.write(fmt::format("{}/init", prefix), spec.init);
    write_uint8(
        writer,
        fmt::format("{}/estimate", prefix),
        static_cast<uint8_t>(spec.estimate));
}

auto write_genetic_spec(
    io::detail::BinaryWriter& writer,
    std::string_view prefix,
    const bayes::GeneticSpec& spec) -> void
{
    write_uint8(
        writer,
        fmt::format("{}/mode", prefix),
        static_cast<uint8_t>(spec.mode));
    write_variance_spec(
        writer, fmt::format("{}/variance", prefix), spec.variance);
    if (spec.sign)
    {
        write_categorical_spec(
            writer, fmt::format("{}/sign", prefix), *spec.sign);
    }
}

auto write_genetic_prior(
    io::detail::BinaryWriter& writer,
    std::string_view prefix,
    const bayes::GeneticPrior& prior) -> void
{
    std::visit(
        [&](const auto& spec)
        {
            using T = std::decay_t<decltype(spec)>;
            if constexpr (std::is_same_v<T, bayes::GeneticSpec>)
            {
                write_uint8(writer, fmt::format("{}/kind", prefix), 0);
                write_genetic_spec(
                    writer, fmt::format("{}/spec", prefix), spec);
            }
            else
            {
                write_uint8(writer, fmt::format("{}/kind", prefix), 1);
                write_genetic_spec(
                    writer,
                    fmt::format("{}/spec/additive", prefix),
                    spec.additive);
                write_genetic_spec(
                    writer,
                    fmt::format("{}/spec/dominance", prefix),
                    spec.dominance);
            }
        },
        prior.spec);

    if (prior.mixture)
    {
        write_uint8(writer, fmt::format("{}/mixture/present", prefix), 1);
        write_uint8(
            writer,
            fmt::format("{}/mixture/strategy", prefix),
            static_cast<uint8_t>(prior.mixture->strategy.index()));
        if (const auto* sm
            = std::get_if<bayes::ScaledMixture>(&prior.mixture->strategy))
        {
            writer.write(
                fmt::format("{}/mixture/multiplier", prefix), sm->multiplier);
        }
        write_categorical_spec(
            writer,
            fmt::format("{}/mixture/proportions", prefix),
            prior.mixture->proportions);
    }
    else
    {
        write_uint8(writer, fmt::format("{}/mixture/present", prefix), 0);
    }
}

auto write_method(
    io::detail::BinaryWriter& writer,
    const bayes::BayesMethod& method) -> void
{
    for (uint8_t i = 0; i < static_cast<uint8_t>(method.genetics.size()); ++i)
    {
        write_genetic_prior(
            writer, fmt::format("method/genetic/{}", i), method.genetics[i]);
    }

    for (uint8_t i = 0; i < static_cast<uint8_t>(method.randoms.size()); ++i)
    {
        write_variance_spec(
            writer, fmt::format("method/random/{}", i), method.randoms[i]);
    }

    write_variance_spec(writer, "method/residual", method.residual);
}

}  // namespace

auto write_checkpoint(
    const mcmc::State& state,
    const std::mt19937_64& rng,
    const bayes::BayesMethod& method,
    std::string_view prefix) -> void
{
    io::detail::BinaryWriter writer(fmt::format("{}.ckpt", prefix));

    write_method(writer, method);
    write_fixed(writer, state.fixed());
    write_random(writer, state.random());
    write_genetics(writer, state.genetics());
    write_residual(writer, state.residual());
    write_rng(writer, rng);
}

}  // namespace gelex
