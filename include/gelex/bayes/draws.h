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

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/detail/draws_factory.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance/draws.h"
#include "gelex/exception.h"
#include "gelex/infra/log.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"

namespace gelex
{

template <typename GeneticPrior>
class BayesDraws;

template <typename GeneticPrior>
[[nodiscard]] auto make_draws(
    const BayesPrior<GeneticPrior>& prior,
    const BayesModel& model,
    std::string_view output_path,
    std::uint64_t draw_count) -> BayesDraws<GeneticPrior>;

class RandomEffectDraws
{
   public:
    RandomEffectDraws(VectorDraw coefficients, ScalarDraw variance)
        : coefficients_{std::move(coefficients)}, variance_{std::move(variance)}
    {
    }

    auto append(const RandomEffectState& state) -> void
    {
        coefficients_.append(state.coefficients);
        variance_.append(state.variance);
    }

    [[nodiscard]] auto coefficients() const noexcept -> const VectorDraw&
    {
        return coefficients_;
    }

    [[nodiscard]] auto variance() const noexcept -> const ScalarDraw&
    {
        return variance_;
    }

   private:
    VectorDraw coefficients_;
    ScalarDraw variance_;
};

template <typename GeneticPrior>
class BayesDraws
{
   public:
    using genetic_draws_type = detail::genetic_draws_t<GeneticPrior>;

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
        genetic_.append(state.genetic());
        residual_.append(state.residual().variance);
        variance_summary_.append(variance_summary);
        ++appended_;
    }

    [[nodiscard]] auto fixed() const noexcept -> const VectorDraw&
    {
        return fixed_;
    }

    [[nodiscard]] auto random() const noexcept
        -> std::span<const RandomEffectDraws>
    {
        return random_;
    }

    [[nodiscard]] auto genetic() const noexcept -> const genetic_draws_type&
    {
        return genetic_;
    }

    [[nodiscard]] auto residual() const noexcept -> const ScalarDraw&
    {
        return residual_;
    }

    [[nodiscard]] auto variance_summary() const noexcept
        -> const VarianceSummaryDraws<GeneticPrior::modes>&
    {
        return variance_summary_;
    }

   private:
    template <typename T>
    friend auto make_draws(
        const BayesPrior<T>& prior,
        const BayesModel& model,
        std::string_view output_path,
        std::uint64_t draw_count) -> BayesDraws<T>;

    BayesDraws(
        const BayesPrior<GeneticPrior>& prior,
        const BayesModel& model,
        std::string_view output_path,
        std::uint64_t draw_count)
        : writer_{output_path},
          fixed_{writer_.reserve<float>(
              "fixed/coefficients",
              BinaryShape{
                  static_cast<std::uint64_t>(model.fixed().X().cols()),
                  draw_count})},
          random_{make_random_draws(model, writer_, draw_count)},
          genetic_{detail::make_draws(
              prior.genetic(),
              model.genetic(),
              writer_,
              draw_count)},
          residual_{writer_.reserve<double>(
              "residual/variance",
              BinaryShape{1, draw_count})},
          variance_summary_{model.random(), writer_, draw_count},
          draw_count_{draw_count}
    {
    }

    static auto make_random_draws(
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
                VectorDraw{writer.reserve<float>(
                    fmt::format("random/{}/coefficients", design.name()),
                    BinaryShape{
                        static_cast<std::uint64_t>(design.X().cols()),
                        draw_count})},
                ScalarDraw{writer.reserve<double>(
                    fmt::format("random/{}/variance", design.name()),
                    BinaryShape{1, draw_count})});
        }
        return random;
    }

    // Every leaf PayloadWriter borrows this address.
    BinaryWriter writer_;
    VectorDraw fixed_;
    std::vector<RandomEffectDraws> random_;
    genetic_draws_type genetic_;
    ScalarDraw residual_;
    VarianceSummaryDraws<GeneticPrior::modes> variance_summary_;
    std::uint64_t draw_count_;
    std::uint64_t appended_{0};
};

template <typename GeneticPrior>
[[nodiscard]] auto make_draws(
    const BayesPrior<GeneticPrior>& prior,
    const BayesModel& model,
    std::string_view output_path,
    std::uint64_t draw_count) -> BayesDraws<GeneticPrior>
{
    return BayesDraws<GeneticPrior>{prior, model, output_path, draw_count};
}

}  // namespace gelex

#endif  // GELEX_BAYES_DRAWS_H_
