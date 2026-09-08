// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/draws.h"

#include <cstdint>
#include <fmt/format.h>
#include <vector>

#include "gelex/bayes/model.h"
#include "gelex/bayes/random_design.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"

namespace gelex::detail
{

auto make_random_draws(
    const BayesModel& model,
    BinaryWriter& writer,
    std::uint64_t draw_count) -> std::vector<RandomEffectDraws>
{
    const auto designs = model.random();
    std::vector<RandomEffectDraws> random;
    random.reserve(designs.size());
    for (const auto& design : designs)
    {
        random.emplace_back(
            writer.reserve<float>(
                fmt::format("random/{}/coefficients", design.name()),
                BinaryShape{
                    static_cast<std::uint64_t>(design.X().cols()), draw_count}),
            writer.reserve<double>(
                fmt::format("random/{}/variance", design.name()),
                BinaryShape{1, draw_count}));
    }
    return random;
}

}  // namespace gelex::detail
