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

#ifndef GELEX_BAYES_DRAWS_H_
#define GELEX_BAYES_DRAWS_H_

#include <cassert>
#include <cstdint>
#include <exception>
#include <fmt/format.h>
#include <ranges>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/construction.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance/draws.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/log.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

class RandomEffectDraws
{
   public:
    RandomEffectDraws(
        PayloadWriter<float> coefficients,
        PayloadWriter<double> variance)
        : coefficients_{std::move(coefficients)}, variance_{std::move(variance)}
    {
    }

    auto append(const RandomEffectState& state) -> void
    {
        coefficients_.append(state.coefficients.cast<float>().eval());
        variance_.append(state.variance);
    }

   private:
    PayloadWriter<float> coefficients_;
    PayloadWriter<double> variance_;
};

GELEX_NAMESPACE_BEGIN(detail)
inline auto make_random_draws(
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
GELEX_NAMESPACE_END(detail)

template <typename GeneticPrior>
class BayesDraws
{
   public:
    using genetic_draws_type = detail::genetic_draws_t<GeneticPrior>;

    BayesDraws(
        const BayesPrior<GeneticPrior>& prior,
        const BayesModel& model,
        std::string_view output_path,
        std::uint64_t draw_count)
        : writer_{output_path},
          fixed_{writer_.reserve<double>(
              "fixed/coefficients",
              BinaryShape{
                  static_cast<std::uint64_t>(model.fixed().X().cols()),
                  draw_count})},
          random_{detail::make_random_draws(model, writer_, draw_count)},
          genetic_{detail::make_draws(
              prior.genetic(),
              model.genetic(),
              writer_,
              draw_count)},
          residual_{writer_.reserve<double>(
              "residual/variance",
              BinaryShape{1, draw_count})},
          variance_summary_{writer_, draw_count},
          draw_count_{draw_count}
    {
    }
    BayesDraws(const BayesDraws&) = delete;
    BayesDraws(BayesDraws&&) = delete;
    auto operator=(const BayesDraws&) -> BayesDraws& = delete;
    auto operator=(BayesDraws&&) -> BayesDraws& = delete;
    ~BayesDraws() noexcept
    {
        if (appended_ == draw_count_ || std::uncaught_exceptions() > 0)
        {
            return;
        }
        try
        {
            warn(
                fmt::format(
                    "recorded {} of {} reserved draws",
                    appended_,
                    draw_count_));
        }
        catch (...)  // NOLINT(bugprone-empty-catch): dtor must be noexcept
        {
        }
    }

    auto append(const BayesState<GeneticPrior>& state) -> void
    {
        assert(random_.size() == state.random().size());
        if (appended_ == draw_count_)
        {
            throw GelexException(
                fmt::format(
                    "draw count exceeded: {} draws reserved", draw_count_));
        }
        const auto variance_summary = make_variance_summary(state);
        fixed_.append(state.fixed().coefficients);
        for (auto&& [draws, random_state] :
             std::views::zip(random_, state.random()))
        {
            draws.append(random_state);
        }
        genetic_.for_each(
            [&]<GeneticMode Mode>(auto& draws)
            { draws.append(state.genetic().template get<Mode>()); });
        if constexpr (requires { genetic_.joint(); })
        {
            genetic_.joint().append(state.genetic().joint());
        }
        residual_.append(state.residual().variance);
        variance_summary_.append(variance_summary);
        ++appended_;
    }

   private:
    // Every leaf PayloadWriter borrows this address.
    BinaryWriter writer_;
    PayloadWriter<double> fixed_;
    std::vector<RandomEffectDraws> random_;
    genetic_draws_type genetic_;
    PayloadWriter<double> residual_;
    VarianceSummaryDraws<GeneticPrior::modes> variance_summary_;
    std::uint64_t draw_count_;
    std::uint64_t appended_{0};
};

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_DRAWS_H_
