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

#include "gelex/model/bayes/writer/checkpoint_writer.h"

#include <fmt/format.h>
#include <cstdint>
#include <span>
#include <sstream>
#include <variant>

#include "gelex/io/binary_writer.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"

namespace gelex
{

namespace
{

auto write_scalar(
    detail::BinaryWriter& writer,
    std::string_view path,
    double value) -> void
{
    auto handle = writer.reserve<double>(path, 1, 1);
    writer.write(handle, value);
}

auto write_uint8(
    detail::BinaryWriter& writer,
    std::string_view path,
    uint8_t value) -> void
{
    auto handle = writer.reserve<uint8_t>(path, 1, 1);
    writer.write(handle, value);
}

auto write_fixed(detail::BinaryWriter& writer, const bayes::FixedState& fs)
    -> void
{
    writer.write(fmt::format("{}/coeff", EffectType::fixed()), fs.coeffs);
}

auto write_random(
    detail::BinaryWriter& writer,
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
    detail::BinaryWriter& writer,
    const EffectType& effect,
    const bayes::Assignment& a) -> void
{
    writer.write(fmt::format("{}/group/assignment", effect), a.tracker);
    writer.write(fmt::format("{}/group/proportion", effect), a.proportion);
    writer.write(fmt::format("{}/group/stick_probs", effect), a.stick_probs);
}

auto write_component_allocation(
    detail::BinaryWriter& writer,
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
    detail::BinaryWriter& writer,
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
    detail::BinaryWriter& writer,
    const EffectType& effect,
    const bayes::Assignment& sign) -> void
{
    writer.write(fmt::format("{}/sign/assignment", effect), sign.tracker);
    writer.write(fmt::format("{}/sign/proportion", effect), sign.proportion);
    writer.write(fmt::format("{}/sign/stick_probs", effect), sign.stick_probs);
}

auto write_genetics(
    detail::BinaryWriter& writer,
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
    detail::BinaryWriter& writer,
    const bayes::ResidualState& rs) -> void
{
    writer.write(fmt::format("{}/adj_pheno", EffectType::residual()), rs.y_adj);
    write_scalar(
        writer,
        fmt::format("{}/variance", EffectType::residual()),
        rs.variance);
}

auto write_rng(detail::BinaryWriter& writer, const std::mt19937_64& rng) -> void
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

auto write_variance_prior(
    detail::BinaryWriter& writer,
    std::string_view prefix,
    const bayes::VariancePrior& vp) -> void
{
    write_scalar(writer, fmt::format("{}/nu", prefix), vp.param.nu);
    write_scalar(writer, fmt::format("{}/s2", prefix), vp.param.s2);
    write_scalar(writer, fmt::format("{}/init", prefix), vp.init);
    write_scalar(
        writer, fmt::format("{}/size", prefix), static_cast<double>(vp.size));
}

auto write_genetic_prior(
    detail::BinaryWriter& writer,
    const bayes::GeneticPrior& gp) -> void
{
    const auto effect = EffectType::from_genetic(gp.type);
    const auto prefix = fmt::format("prior/{}", effect);

    const auto type_index = static_cast<uint8_t>(gp.marker.index());
    write_uint8(writer, fmt::format("{}/type", prefix), type_index);

    std::visit(
        [&](const auto& mp)
        {
            write_variance_prior(
                writer, fmt::format("{}/variance", prefix), mp.variance);
            using T = std::decay_t<decltype(mp)>;
            if constexpr (
                std::is_same_v<T, bayes::SpikePrior>
                || std::is_same_v<T, bayes::MixturePrior>)
            {
                writer.write(
                    fmt::format("{}/proportion/init", prefix),
                    mp.proportion.init);
                write_uint8(
                    writer,
                    fmt::format("{}/proportion/estimate", prefix),
                    static_cast<uint8_t>(mp.proportion.estimate));
            }
            if constexpr (std::is_same_v<T, bayes::MixturePrior>)
            {
                writer.write(
                    fmt::format("{}/multiplier", prefix), mp.multiplier);
            }
        },
        gp.marker);

    if (gp.sign)
    {
        write_scalar(
            writer, fmt::format("{}/sign", prefix), gp.sign->init_value);
    }
}

auto write_priors(detail::BinaryWriter& writer, const bayes::Priors& priors)
    -> void
{
    write_variance_prior(
        writer,
        fmt::format("prior/{}", EffectType::residual()),
        priors.residual());

    for (size_t i = 0; i < priors.random().size(); ++i)
    {
        write_variance_prior(
            writer,
            fmt::format("prior/{}/{}", EffectType::random(), i),
            priors.random()[i]);
    }

    for (const auto& gp : priors.genetics())
    {
        write_genetic_prior(writer, gp);
    }
}

}  // namespace

auto write_checkpoint(
    const mcmc::State& state,
    const std::mt19937_64& rng,
    const bayes::Priors& priors,
    std::string_view prefix) -> void
{
    detail::BinaryWriter writer(fmt::format("{}.ckpt", prefix));

    write_priors(writer, priors);
    write_fixed(writer, state.fixed());
    write_random(writer, state.random());
    write_genetics(writer, state.genetics());
    write_residual(writer, state.residual());
    write_rng(writer, rng);
}

}  // namespace gelex
